#!/usr/bin/env nextflow

/*
 *   Nextclade consensus FASTA workflow
 *   ------------------------------------
 *   Takes pre-assembled consensus FASTAs through QC, per-segment Nextclade alignment
 *   against the --db bundle's datasets, and conversion of the JSON to a freebayes-style VCF.
 *
 *   The sample name comes from the input filename (modules/shared/sample_id.nf); input
 *   headers supply only the segment label. Accepted header formats:
 *     Underscore-delimited (legacy): >SAMPLE_ID_segmentN_SEGMENT
 *       e.g. >WIN-AH-2024-FAV-0415-1-OS_segment1_PB2
 *     Pipe-delimited (GISAID-style): >EPI_ISL_ACCESSION|subtype|SEGMENT|...|clade_SEGMENT
 *       e.g. >EPI_ISL_20399898|A_/_H5N1|NA|6|EPI5261179|2.3.4.4b_NA
 *   Cleaned headers are normalised to {sample_id}_{SEGMENT} by prepare_consensus.
 */

include { hash_files }                           from '../modules/shared/hash_files.nf'
include { prepare_consensus }                    from '../modules/shared/consensus_qc.nf'
include { flumut }                               from '../modules/shared/flumut.nf'
include { nextclade_run_batch }                  from '../modules/nextclade/nextclade.nf'
include { coordinate_conversion }                from '../modules/nextclade/coordinate_conversion.nf'
include { merge_flumut_mutations }               from '../modules/nextclade/merge_flumut_mutations.nf'
include { mutation_watch }                       from '../modules/shared/mutation_watch.nf'
include { pipeline_provenance }                  from '../modules/shared/provenance.nf'
include { collect_provenance }                   from '../modules/shared/provenance.nf'
include { resolve_db }                           from '../modules/shared/database.nf'
include { consensus_fasta_channel }              from '../modules/shared/sample_id.nf'
include { samplesheet_channel }                  from '../modules/shared/sample_id.nf'

workflow NEXTCLADE {

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

    // Consensus FASTAs from the sample sheet or the input directory.
    // Each file is expected to contain 8 segments for a single sample.
    ch_fasta_input = params.samplesheet_input != 'NO_FILE'
        ? samplesheet_channel(params.samplesheet_input, ['FASTA'])
        : consensus_fasta_channel(params.fasta_input)

    main:

        // Hash the input FASTA. One file per sample, so no list-wrapping here.
        hash_files(ch_fasta_input.combine(Channel.of("fasta_input")))

        // --- Step 1: Validate and clean input consensus FASTA ---
        prepare_consensus(ch_fasta_input)

        // One run-level table of the per-segment consensus QC rows. The name carries the
        // workflow because main.nf runs both FASTA branches into the same outdir, and two
        // collectFile operations must not target one path.
        prepare_consensus.out.qc
            .map{ it -> it[1] }
            .collectFile(
                name: 'nextclade_consensus_qc_summary.tsv',
                storeDir: "${params.outdir}/aggregate",
                keepHeader: true,                   // keep header from first file only
                skip: 1,                            // drop the header line from the rest
                sort: false
            )

        // --- Step 2: HPAI mutation screening (per sample, parallel branch) ---
        flumut(prepare_consensus.out.cleaned)

        // --- Step 3: Collect each segment across all samples into one FASTA ---
        //   Cleaned headers are {sample_id}_{SEGMENT} with a canonical SEGMENT; the batch
        //   header is the bare sample_id, which nextclade reports back as seqName.
        ch_segment_batches = prepare_consensus.out.cleaned
            .splitFasta(elem: 1, record: [id: true, seqString: true])     // [sample_id, record]
            .collectFile { item ->
                def (sample_id, record) = item
                ["${record.id - "${sample_id}_"}.fasta", ">${sample_id}\n${record.seqString}\n"]
            }
            .map { fasta -> [fasta.baseName, fasta, db.nextclade.resolve(fasta.baseName)] }
            // -> [segment, batch.fasta, dataset_dir]

        // --- Step 4: Run nextclade once per segment for all samples, and convert to VCF ---
        //   Also writes each sample's report-only consensus/reference identity QC.
        //   With --emit_no_mature_vcf it also publishes each VCF without mature-peptide
        //   (HA1/HA2) annotations. Nothing downstream reads that copy: mutation_watch and
        //   merge_flumut_mutations resolve HA markers through the HA1/HA2 Feature_ID.
        nextclade_run_batch(ch_segment_batches)

        // --- Step 5: Scan the final VCFs for custom watchlists (optional) ---
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

        // --- Steps 6-7: Coordinate conversion + flumut mutation merging (optional) ---
        //   Gated behind --merge_flumut_mutations (default false); otherwise the workflow
        //   stops at the *_final.vcf produced above.
        if (params.merge_flumut_mutations) {

            // --- Step 6: Build the flumut-to-Nextclade coordinate map ---
            //   source = the flumutdb package's reference (lift FROM), target = the --db
            //   datasets (lift ONTO). The target label comes from the bundle's db.yml, so the
            //   published map names the reference rather than where the bundle sits.
            coordinate_conversion(Channel.value([db.nextclade_ref_name, db.nextclade]))

            // --- Step 7: Annotate every sample's VCFs with FluMutDB markers, one task per segment ---
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
        ch_provenance = ch_fasta_input.map{ it -> it[0] }
        ch_provenance = ch_provenance.combine(ch_pipeline_provenance).map{ it -> [it[0], [it[1]]] }
        ch_provenance = ch_provenance.join(hash_files.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
        ch_provenance = ch_provenance.join(prepare_consensus.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
        ch_provenance = ch_provenance.join(flumut.out.provenance).map{ it -> [it[0], it[1] << it[2]] }

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
