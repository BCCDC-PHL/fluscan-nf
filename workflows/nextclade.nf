#!/usr/bin/env nextflow

/*
 *   Nextclade consensus FASTA workflow
 *   ------------------------------------
 *   Takes the consensus FASTAs FASTA_PREP has cleaned through per-segment Nextclade
 *   alignment against each sample's bundle's datasets, and conversion of the JSON to a
 *   freebayes-style VCF.
 */

include { nextclade_run_batch }                  from '../modules/nextclade/nextclade.nf'
include { coordinate_conversion }                from '../modules/nextclade/coordinate_conversion.nf'
include { merge_flumut_mutations }               from '../modules/nextclade/merge_flumut_mutations.nf'
include { mutation_watch }                       from '../modules/shared/mutation_watch.nf'
include { pipeline_provenance }                  from '../modules/shared/provenance.nf'
include { collect_provenance }                   from '../modules/shared/provenance.nf'
include { resolve_dbs }                          from '../modules/shared/database.nf'
include { input_db_paths }                       from '../modules/shared/sample_id.nf'

workflow NEXTCLADE {

    take:
        ch_cleaned              // [sample_id, cleaned consensus FASTA]
        ch_db                   // [sample_id, bundle name]
        ch_prep_provenance      // [sample_id, [hash, prepare, flumut]]

    main:

        // Only the Nextclade datasets are read here; this branch never BLASTs or runs SnpEff.
        def dbs = resolve_dbs(input_db_paths(), ['nextclade'])

        ch_workflow_metadata = Channel.fromList(dbs.values().collect { db -> [
            db.name,
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
        ] })

        ch_pipeline_provenance = pipeline_provenance(ch_workflow_metadata)    // [bundle name, record]

        // --- Step 1: Collect each segment across one bundle's samples into one FASTA ---
        //   Cleaned headers are {sample_id}_{SEGMENT} with a canonical SEGMENT; the batch
        //   header is the bare sample_id, which nextclade reports back as seqName.
        ch_segment_batches = ch_cleaned.join(ch_db)
            .splitFasta(elem: 1, record: [id: true, seqString: true])     // [sample_id, record, bundle name]
            .collectFile { item ->
                def (sample_id, record, db_name) = item
                ["${db_name}__${record.id - "${sample_id}_"}.fasta", ">${sample_id}\n${record.seqString}\n"]
            }
            .map { fasta ->
                def (db_name, segment) = (fasta.baseName =~ /^(.+)__([^_]+)$/)[0][1..2]
                [db_name, segment, fasta, dbs[db_name].nextclade.resolve(segment)]
            }
            // -> [bundle name, segment, batch.fasta, dataset_dir]

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
            // [bundle, segment, [vcfs]] -> [sample_id, [segments], [vcfs]]: one scan task per sample.
            ch_sample_vcfs = nextclade_run_batch.out.vcf
                .transpose()
                .map { db_name, segment, vcf -> [vcf.name.split('__')[0], segment, vcf] }
                .groupTuple()
            mutation_watch(ch_sample_vcfs, ch_watchlists.collect())
        }

        // --- Steps 4-5: Coordinate conversion + flumut mutation merging (optional) ---
        //   Gated behind --merge_flumut_mutations (default false); otherwise the workflow
        //   stops at the *_final.vcf produced above.
        if (params.merge_flumut_mutations) {

            // --- Step 4: Build the flumut-to-Nextclade coordinate map, once per bundle ---
            //   source = the flumutdb package's reference (lift FROM), target = the bundle's
            //   datasets (lift ONTO). The target label comes from the bundle's db.yml, so the
            //   published map names the reference rather than where the bundle sits.
            coordinate_conversion(Channel.fromList(dbs.values().collect { db -> [db.name, db.nextclade_ref_name, db.nextclade] }))

            // --- Step 5: Annotate every sample's VCFs with FluMutDB markers, one task per segment ---
            //   Additive: writes a separate *_flumut_annotated.vcf and leaves *_final.vcf
            //   untouched. Each task stages only its own segment's shard.

            // Key each per-segment shard by its bundle and the segment label parsed from the filename.
            ch_shards = coordinate_conversion.out.json_shards
                .transpose()
                .map { db_name, f ->                           // [bundle, SEG, shard]
                    def matcher = f.name =~ /coordinate_map\.([^.]+)\.json/
                    if (!matcher) {
                        error "coordinate_conversion emitted '${f.name}', which does not follow " +
                              "coordinate_map.<SEGMENT>.json; cannot key it by segment"
                    }
                    [db_name, matcher[0][1], f]
                }

            // [bundle, segment, [vcfs]] + shard -> [bundle, segment, [vcfs], coord_shard]
            merge_flumut_mutations(nextclade_run_batch.out.vcf.combine(ch_shards, by: [0, 1]))
        }

        // Per-sample provenance, keyed [bundle, sample_id, records] until the end. Nextclade
        // runs once per bundle and segment, so each segment's record is combined onto every
        // sample of that bundle.
        ch_nextclade_provenance = nextclade_run_batch.out.provenance
            .map { db_name, segment, record -> [db_name, record] }
            .groupTuple()
        ch_provenance = ch_prep_provenance.join(ch_db)
            .map { sample_id, records, db_name -> [db_name, sample_id, records] }
            .combine(ch_pipeline_provenance, by: 0)
            .combine(ch_nextclade_provenance, by: 0)
            .map { db_name, sample_id, records, pipeline, nextclade -> [db_name, sample_id, [pipeline] + records + nextclade] }

        if (params.watchlist_dir != 'NO_FILE') {
            // One scan task per sample. remainder: true keeps samples whose scan failed.
            ch_provenance = ch_provenance
                .map { db_name, sample_id, records -> [sample_id, records] }
                .join(mutation_watch.out.provenance, remainder: true)
                .map { sample_id, records, watch -> [sample_id, (records ?: []) + (watch ? [watch] : [])] }
                .join(ch_db)
                .map { sample_id, records, db_name -> [db_name, sample_id, records] }
        }

        if (params.merge_flumut_mutations) {
            // One merge task per bundle and segment, so each segment's record is combined onto
            // every sample of that bundle. A bundle with no VCF to merge gets an empty list.
            ch_merge_provenance = merge_flumut_mutations.out.provenance
                .map { db_name, segment, record -> [db_name, record] }
                .groupTuple()
            ch_merge_provenance = Channel.fromList(dbs.keySet().collect { [it] })
                .join(ch_merge_provenance, remainder: true)
                .map { db_name, records -> [db_name, records ?: []] }

            // coordinate_conversion runs once per bundle, not per sample, so its provenance
            // file is combined onto every sample of that bundle.
            ch_provenance = ch_provenance
                .combine(ch_merge_provenance, by: 0)
                .combine(coordinate_conversion.out.provenance, by: 0)
                .map { db_name, sample_id, records, merge, coord -> [db_name, sample_id, records + merge + [coord]] }
        }

        collect_provenance(ch_provenance.map { db_name, sample_id, records -> [sample_id, records] })

}
