#!/usr/bin/env nextflow

/*
 *   Nextclade consensus FASTA workflow
 *   ------------------------------------
 *   Takes the consensus FASTAs FASTA_PREP has cleaned through per-segment Nextclade
 *   alignment against the --db bundle's datasets, and conversion of the JSON to a
 *   freebayes-style VCF.
 */

include { nextclade_run_batch }                  from '../modules/nextclade/nextclade.nf'
include { coordinate_conversion }                from '../modules/nextclade/coordinate_conversion.nf'
include { merge_flumut_mutations }               from '../modules/nextclade/merge_flumut_mutations.nf'
include { mutation_watch }                       from '../modules/shared/mutation_watch.nf'
include { pipeline_provenance }                  from '../modules/shared/provenance.nf'
include { collect_provenance }                   from '../modules/shared/provenance.nf'
include { resolve_db }                           from '../modules/shared/database.nf'

workflow NEXTCLADE {

    take:
        ch_cleaned              // [sample_id, cleaned consensus FASTA]
        ch_prep_provenance      // [sample_id, [hash, prepare, flumut]]

    main:

        // Only the Nextclade datasets are read here; this branch never BLASTs or runs SnpEff.
        def db = resolve_db(params.db, ['nextclade'])

        ch_workflow_metadata = Channel.value([
            workflow.sessionId,
            workflow.runName,
            workflow.manifest.name,
            workflow.manifest.version,
            workflow.start,
            [
                database_name:    db.name,
                database_version: db.version,
                database_path:    db.root,
            ],
        ])

        ch_pipeline_provenance = pipeline_provenance(ch_workflow_metadata)

        // --- Step 1: Collect each segment across all samples into one FASTA ---
        //   Cleaned headers are {sample_id}_{SEGMENT} with a canonical SEGMENT; the batch
        //   header is the bare sample_id, which nextclade reports back as seqName.
        ch_segment_batches = ch_cleaned
            .splitFasta(elem: 1, record: [id: true, seqString: true])     // [sample_id, record]
            .collectFile { item ->
                def (sample_id, record) = item
                ["${record.id - "${sample_id}_"}.fasta", ">${sample_id}\n${record.seqString}\n"]
            }
            .map { fasta -> [fasta.baseName, fasta, db.nextclade.resolve(fasta.baseName)] }
            // -> [segment, batch.fasta, dataset_dir]

        // --- Step 2: Run nextclade once per segment for all samples, and convert to VCF ---
        //   Also writes each sample's report-only consensus/reference identity QC.
        //   With --emit_no_mature_vcf it also publishes each VCF without mature-peptide
        //   (HA1/HA2) annotations. Nothing downstream reads that copy: mutation_watch and
        //   merge_flumut_mutations resolve HA markers through the HA1/HA2 Feature_ID.
        nextclade_run_batch(ch_segment_batches)

        // --- Step 3: Scan the final VCFs for custom watchlists (optional) ---
        //   Independent of --merge_flumut_mutations: FLUMUT= is a curated annotation written
        //   into the VCF, this is a custom list scanned out of it into a separate TSV.
        if (params.watchlist_dir != 'NO_FILE') {
            // checkIfExists so a mistyped --watchlist_dir fails loudly rather than
            // silently never scanning.
            ch_watchlists = Channel.fromPath(
                params.watchlist_dir + "/*watchlist.csv", checkIfExists: true
            )
            // [segment, [vcfs]] -> [sample_id, [segments], [vcfs]]: one scan task per sample.
            ch_sample_vcfs = nextclade_run_batch.out.vcf
                .transpose()
                .map { segment, vcf -> [vcf.name.split('__')[0], segment, vcf] }
                .groupTuple()
            mutation_watch(ch_sample_vcfs, ch_watchlists.collect())
        }

        // --- Steps 4-5: Coordinate conversion + flumut mutation merging (optional) ---
        //   Gated behind --merge_flumut_mutations (default false); otherwise the workflow
        //   stops at the *_final.vcf produced above.
        if (params.merge_flumut_mutations) {

            // --- Step 4: Build the flumut-to-Nextclade coordinate map ---
            //   source = the flumutdb package's reference (lift FROM), target = the --db
            //   datasets (lift ONTO). The target label comes from the bundle's db.yml, so the
            //   published map names the reference rather than where the bundle sits.
            coordinate_conversion(Channel.value([db.nextclade_ref_name, db.nextclade]))

            // --- Step 5: Annotate every sample's VCFs with FluMutDB markers, one task per segment ---
            //   Additive: writes a separate *_flumut_annotated.vcf and leaves *_final.vcf
            //   untouched. Each task stages only its own segment's shard.

            // Key each per-segment shard by its segment label parsed from the filename.
            ch_shards = coordinate_conversion.out.json_shards
                .flatten()
                .map { f ->                                    // [SEG, shard]
                    def matcher = f.name =~ /coordinate_map\.([^.]+)\.json/
                    if (!matcher) {
                        error "coordinate_conversion emitted '${f.name}', which does not follow " +
                              "coordinate_map.<SEGMENT>.json; cannot key it by segment"
                    }
                    [matcher[0][1], f]
                }

            // [segment, [vcfs]] + shard -> [segment, [vcfs], coord_shard]
            merge_flumut_mutations(nextclade_run_batch.out.vcf.combine(ch_shards, by: 0))
        }

        // Per-sample provenance. The Nextclade steps are keyed by segment rather than by
        // sample, so they have no per-sample record to join here.
        ch_provenance = ch_prep_provenance
            .combine(ch_pipeline_provenance)
            .map{ sample_id, records, pipeline -> [sample_id, [pipeline] + records] }

        if (params.watchlist_dir != 'NO_FILE') {
            // One scan task per sample. remainder: true keeps samples whose scan failed.
            ch_provenance = ch_provenance
                .join(mutation_watch.out.provenance, remainder: true)
                .map{ sample_id, records, watch -> [sample_id, (records ?: []) + (watch ? [watch] : [])] }
        }

        if (params.merge_flumut_mutations) {
            // One merge task per segment, so every segment's record is combined onto every
            // sample. ifEmpty keeps the chain alive when no segment produced a VCF.
            ch_merge_provenance = merge_flumut_mutations.out.provenance
                .map { segment, record -> record }
                .collect()
                .ifEmpty([])
                .map { records -> [records] }
            ch_provenance = ch_provenance
                .combine(ch_merge_provenance)
                .map{ sample_id, records, merge -> [sample_id, records + merge] }

            // coordinate_conversion runs once per pipeline invocation, not per sample, so
            // its single provenance file is combined onto every sample rather than joined.
            ch_provenance = ch_provenance
                .combine(coordinate_conversion.out.provenance)
                .map{ it -> [it[0], it[1] << it[2]] }
        }

        collect_provenance(ch_provenance)

}
