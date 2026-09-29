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
include { split_fasta as split_fasta_consensus } from '../modules/shared/consensus_qc.nf'
include { consensus_identity_qc }                from '../modules/nextclade/consensus_qc.nf'
include { flumut }                               from '../modules/shared/flumut.nf'
include { nextclade_run_batch }                  from '../modules/nextclade/nextclade.nf'
include { nextclade_to_vcf_batch }               from '../modules/nextclade/nextclade.nf'
include { post_process_vcf_batch }               from '../modules/shared/post_process.nf'
include { post_process_vcf_batch as post_process_vcf_batch_no_mature } from '../modules/shared/post_process.nf'
include { coordinate_conversion }                from '../modules/nextclade/coordinate_conversion.nf'
include { merge_flumut_mutations }               from '../modules/nextclade/merge_flumut_mutations.nf'
include { mutation_watch }                       from '../modules/shared/mutation_watch.nf'
include { pipeline_provenance }                  from '../modules/shared/provenance.nf'
include { collect_provenance }                   from '../modules/shared/provenance.nf'
include { resolve_db }                           from '../modules/shared/database.nf'
include { db_segments }                          from '../modules/shared/database.nf'
include { consensus_fasta_channel }              from '../modules/shared/sample_id.nf'
include { samplesheet_channel }                  from '../modules/shared/sample_id.nf'

workflow NEXTCLADE {

    // The canonical segment set, and the only labels coordinate_conversion shards the map by.
    def SEGMENTS = db_segments()

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

    // A name that does not match yields 'unknown_segment' rather than throwing, so the
    // unknown-segment branch below can report it.
    def extract_segment = { sample_id, fasta ->
        def matcher = fasta.name =~ /.*_([A-Za-z0-9-]+)\.s.*/
        def segment = matcher ? matcher[0][1] : "unknown_segment"
        return [sample_id, segment, fasta]
    }

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

        // Split once; reused by the identity QC and the nextclade batching below.
        ch_consensus_segments = split_fasta_consensus(prepare_consensus.out.cleaned)
            .transpose()
            .map(extract_segment)                          // [sample_id, segment, fa]

        // --- Step 1b: Report-only consensus/reference identity QC ---
        ch_identity_qc = ch_consensus_segments
            .map { sample_id, segment, consensus ->
                def ref_fasta = db.nextclade.resolve("${segment}/reference.fasta")
                return [sample_id, segment, consensus, ref_fasta]
            }
        consensus_identity_qc(ch_identity_qc)

        // --- Step 2: HPAI mutation screening (per sample, parallel branch) ---
        flumut(prepare_consensus.out.cleaned)

        // --- Step 3: Group per-segment FASTAs across all samples for batching ---
        ch_segment_batches = ch_consensus_segments
            .map { sample_id, segment, fa -> [segment, sample_id, fa] }
            .groupTuple()                                  // [segment, [sample_ids], [fastas]]
            .map { segment, sample_ids, fastas ->
                def dataset_dir = db.nextclade.resolve(segment)
                return [segment, sample_ids, fastas, dataset_dir]
            }

        // --- Step 4: Run nextclade once per segment for all samples ---
        nextclade_run_batch(ch_segment_batches)

        // --- Step 5: Convert per-sample nextclade JSON to freebayes-style VCF ---
        nextclade_to_vcf_batch(nextclade_run_batch.out.json_with_ref)

        // --- Step 6: Post-process VCFs (fix MNV codons); one task per segment ---
        post_process_vcf_batch(nextclade_to_vcf_batch.out.vcf)

        // [segment, [vcfs]] -> [sample_id, segment, vcf]. Shared by the watchlist scan below
        // and, when enabled, by the flumut merge further down.
        ch_final_vcf = post_process_vcf_batch.out.vcf
            .transpose()
            .map { segment, vcf ->
                def sample_id = vcf.name.split('__')[0]
                [sample_id, segment, vcf]
            }

        // --- Step 6a: The same VCFs without mature-peptide (HA1/HA2) annotations (optional) ---
        //   A parallel output for consumers that require HA annotated as HA0 only. Nothing
        //   downstream in this workflow reads it: mutation_watch and merge_flumut_mutations
        //   resolve HA markers through the HA1/HA2 Feature_ID and stay on ch_final_vcf.
        if (params.emit_no_mature_vcf) {
            post_process_vcf_batch_no_mature(nextclade_to_vcf_batch.out.vcf_no_mature)

            ch_final_vcf_no_mature = post_process_vcf_batch_no_mature.out.vcf
                .transpose()
                .map { segment, vcf ->
                    def sample_id = vcf.name.split('__')[0]
                    [sample_id, segment, vcf]
                }
        }

        // --- Step 6b: Scan the final VCF for a custom watchlist (optional) ---
        //   Independent of --merge_flumut_mutations: FLUMUT= is a curated annotation written
        //   into the VCF, this is a custom list scanned out of it into a separate TSV.
        if (params.watchlist_dir != 'NO_FILE') {
            // checkIfExists so a mistyped --watchlist_dir fails loudly rather than
            // silently never scanning.
            ch_watchlists = Channel.fromPath(
                params.watchlist_dir + "/*watchlist.csv", checkIfExists: true
            )
            mutation_watch(ch_final_vcf.combine(ch_watchlists))
        }

        // --- Steps 7-8: Coordinate conversion + flumut mutation merging (optional) ---
        //   Gated behind --merge_flumut_mutations (default false); otherwise the workflow
        //   stops at the *_final.vcf produced above.
        if (params.merge_flumut_mutations) {

            // --- Step 7: Build the flumut-to-Nextclade coordinate map ---
            //   source = the flumutdb package's reference (lift FROM), target = the --db
            //   datasets (lift ONTO). The target label comes from the bundle's db.yml, so the
            //   published map names the reference rather than where the bundle sits.
            coordinate_conversion(Channel.value([db.nextclade_ref_name, db.nextclade]))

            // --- Step 8: Annotate each per-sample/per-segment VCF with FluMutDB markers ---
            //   Additive: writes a separate *_flumut_annotated.vcf and leaves *_final.vcf
            //   untouched. Run per (sample, segment) so each task stages only its own shard.
            ch_vcf_flat = ch_final_vcf                         // [sample_id, segment, vcf]

            // The combine below matches on the segment label and ch_shards only carries the
            // eight canonical ones. Split first so a VCF labelled anything else is reported
            // rather than vanishing with no merged counterpart and no explanation.
            ch_vcf_by_segment = ch_vcf_flat.branch { sample_id, segment, vcf ->
                known:   SEGMENTS.contains(segment)
                unknown: true
            }

            ch_vcf_by_segment.unknown.subscribe { sample_id, segment, vcf ->
                log.warn "merge_flumut_mutations: ${sample_id} has no coordinate map for segment " +
                         "'${segment}' (not one of ${SEGMENTS.join(', ')}); ${vcf.name} is left " +
                         "unannotated and produces no *_flumut_annotated.vcf"
            }

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

            ch_merge_in = ch_vcf_by_segment.known
                .map { sample_id, segment, vcf -> [segment, sample_id, vcf] }
                .combine(ch_shards, by: 0)                     // + coord shard (by segment)
                .map { segment, sample_id, vcf, shard -> [sample_id, segment, vcf, shard] }
                // -> [sample_id, segment, vcf, coord_shard]

            merge_flumut_mutations(ch_merge_in)
        }

        // Per-sample provenance. post_process_vcf_batch is keyed by segment rather than by
        // sample, so it has no per-sample record to join here.
        ch_provenance = ch_fasta_input.map{ it -> it[0] }
        ch_provenance = ch_provenance.combine(ch_pipeline_provenance).map{ it -> [it[0], [it[1]]] }
        ch_provenance = ch_provenance.join(hash_files.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
        ch_provenance = ch_provenance.join(prepare_consensus.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
        ch_provenance = ch_provenance.join(flumut.out.provenance).map{ it -> [it[0], it[1] << it[2]] }

        if (params.watchlist_dir != 'NO_FILE') {
            // Group one task per (sample, segment, watchlist) back to one entry per sample.
            // remainder: true keeps samples whose scans all failed.
            ch_provenance = ch_provenance
                .join(mutation_watch.out.provenance.groupTuple(), remainder: true)
                .map{ it -> [it[0], (it[1] ?: []) + (it[2] ?: [])] }
        }

        if (params.merge_flumut_mutations) {
            // Group one merge task per (sample, segment) back to one entry per sample.
            // remainder: true keeps samples whose merges all failed.
            ch_provenance = ch_provenance
                .join(merge_flumut_mutations.out.provenance.groupTuple(), remainder: true)
                .map{ it -> [it[0], (it[1] ?: []) + (it[2] ?: [])] }

            // coordinate_conversion runs once per pipeline invocation, not per sample, so
            // its single provenance file is combined onto every sample rather than joined.
            ch_provenance = ch_provenance
                .combine(coordinate_conversion.out.provenance)
                .map{ it -> [it[0], it[1] << it[2]] }
        }

        collect_provenance(ch_provenance)

}
