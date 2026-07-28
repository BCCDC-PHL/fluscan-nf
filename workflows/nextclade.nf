#!/usr/bin/env nextflow

/*
 *   Nextclade consensus FASTA workflow
 *   ------------------------------------
 *   Processes pre-assembled consensus FASTA files through quality control,
 *   per-segment Nextclade alignment against GooseGuangdong1996 reference
 *   datasets, and conversion of Nextclade JSON output to a freebayes-style VCF.
 *
 *   Accepted input FASTA header formats:
 *     Underscore-delimited (legacy): >SAMPLE_ID_segmentN_SEGMENT
 *       e.g. >WIN-AH-2024-FAV-0415-1-OS_segment1_PB2
 *     Pipe-delimited (GISAID-style): >EPI_ISL_ACCESSION|subtype|SEGMENT|...|clade_SEGMENT
 *       e.g. >EPI_ISL_20399898|A_/_H5N1|NA|6|EPI5261179|2.3.4.4b_NA
 *   Cleaned headers are normalised to {sample_id}_{SEGMENT} by prepare_consensus.
 *
 *   Segment datasets are resolved from --nextclade_datasets_dir by lowercasing
 *   the segment label (e.g. HA -> ha, M -> m).
 *
 *   Batching strategy: sequences are grouped by segment across all samples so
 *   nextclade runs 8 times total (once per segment) rather than 8 × N times.
 *   Per-sample JSON/TSV/VCF outputs are produced by splitting the batch results.
 */

nextflow.enable.dsl = 2

include { prepare_consensus }                    from '../modules/shared/consensus_qc.nf'
include { split_fasta as split_fasta_consensus } from '../modules/shared/consensus_qc.nf'
include { consensus_identity_qc }                from '../modules/nextclade/consensus_qc.nf'
include { flumut }                               from '../modules/shared/flumut.nf'
include { nextclade_run_batch }                  from '../modules/nextclade/nextclade.nf'
include { nextclade_to_vcf_batch }               from '../modules/nextclade/nextclade.nf'
include { post_process_vcf_batch }               from '../modules/shared/post_process.nf'
include { coordinate_conversion }                from '../modules/nextclade/coordinate_conversion.nf'
include { merge_flumut_mutations }               from '../modules/nextclade/merge_flumut_mutations.nf'
include { pipeline_provenance }                  from '../modules/shared/provenance.nf'
include { collect_provenance }                   from '../modules/shared/provenance.nf'

workflow NEXTCLADE {

    // The canonical segment set, and the only labels coordinate_conversion shards the map by.
    def SEGMENTS = ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'M', 'NS']

    ch_workflow_metadata = Channel.value([
        workflow.sessionId,
        workflow.runName,
        workflow.manifest.name,
        workflow.manifest.version,
        workflow.start,
    ])

    ch_pipeline_provenance = pipeline_provenance(ch_workflow_metadata)

    def extract_segment = { sample_id, fasta ->
        def segment = fasta =~ /.*_([A-z0-9-]+)\.s.*/
        segment = segment[0][1] ? segment[0][1] : "unknown_segment"
        return [sample_id, segment, fasta]
    }

    // Discover consensus FASTA files in the input directory.
    // Each file is expected to contain 8 segments for a single sample.
    ch_fasta_input = Channel.fromPath(
        params.fasta_input + "/*.{fasta,fa,fna}",
        checkIfExists: true
    ).map { fasta ->
        def sample_id = fasta.baseName.replaceAll(/[_\.]consensus$/, '').replaceAll(/_/,'-')
        sample_id = sample_id.replaceFirst(/^\d{2}-\d+(?:-IP\d+)?-/, "")
        return [sample_id, fasta]
    }

    main:

        // Provenance channel starts with just the sample IDs; per-process YAMLs are joined
        // on as they are produced, and concatenated at the end by collect_provenance.
        ch_provenance = ch_fasta_input.map{ it -> it[0] }

        if (
            !params.nextclade_datasets_dir ||
            params.nextclade_datasets_dir == 'NO_FILE' ||
            !file(params.nextclade_datasets_dir).exists()
        ) {
            error "ERROR: --nextclade_datasets_dir must be set to an existing path for the Nextclade workflow."
        }

        if (
            params.merge_flumut_mutations &&
            (
                !params.flumut_reference_dir ||
                params.flumut_reference_dir == 'NO_FILE' ||
                !file(params.flumut_reference_dir).exists()
            )
        ) {
            error "ERROR: --flumut_reference_dir must be set to an existing path (flumut's Nextclade dataset directory, the source coordinate reference) when --merge_flumut_mutations is enabled."
        }

        // --- Step 1: Validate and clean input consensus FASTA ---
        prepare_consensus(ch_fasta_input)

        // Split each cleaned consensus into per-segment FASTAs once; reused by
        // the identity QC and the cross-sample nextclade batching below.
        ch_consensus_segments = split_fasta_consensus(prepare_consensus.out.cleaned)
            .transpose()
            .map(extract_segment)                          // [sample_id, segment, fa]

        // --- Step 1b: Report-only consensus/reference identity QC ---
        //   Pair each consensus segment with its Nextclade dataset reference.
        ch_identity_qc = ch_consensus_segments
            .map { sample_id, segment, consensus ->
                def ref_fasta = file("${params.nextclade_datasets_dir}/${segment.toLowerCase()}/reference.fasta")
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
                def dataset_dir = file("${params.nextclade_datasets_dir}/${segment.toLowerCase()}")
                return [segment, sample_ids, fastas, dataset_dir]
            }

        // --- Step 4: Run nextclade once per segment for all samples ---
        nextclade_run_batch(ch_segment_batches)

        // --- Step 5: Convert per-sample nextclade JSON to freebayes-style VCF ---
        nextclade_to_vcf_batch(nextclade_run_batch.out.json_with_ref)

        // --- Step 6: Post-process VCFs (fix MNV codons); one task per segment ---
        post_process_vcf_batch(nextclade_to_vcf_batch.out.vcf)

        // Per-sample provenance. post_process_vcf_batch is keyed by segment rather than by
        // sample, so it has no per-sample record to join here.
        ch_provenance = ch_provenance.combine(ch_pipeline_provenance).map{ it -> [it[0], [it[1]]] }
        ch_provenance = ch_provenance.join(prepare_consensus.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
        ch_provenance = ch_provenance.join(flumut.out.provenance).map{ it -> [it[0], it[1] << it[2]] }

        // --- Steps 7-8: Coordinate conversion + flumut mutation merging (optional) ---
        //   Gated behind --merge_flumut_mutations (default false). When enabled,
        //   --flumut_reference_dir is validated above and these two steps run;
        //   otherwise the workflow stops at the *_final.vcf produced above.
        if (params.merge_flumut_mutations) {

            // --- Step 7: Build the flumut-to-Nextclade coordinate map ---
            //   Runs once per pipeline invocation (independent of samples). Translates
            //   nucleotide/amino-acid indices from flumut's fixed reference (source)
            //   onto the selected Nextclade dataset (target) coordinate system, so a
            //   downstream step can remap flumut mutations into Nextclade coordinates.
            //   Idempotent, including the trivial identity case where the selected
            //   Nextclade dataset is derived from flumut's reference.
            //   Note: source = flumut (lift FROM), target = nextclade (lift ONTO).
            ch_source_reference = Channel.value(
                [
                    file(params.flumut_reference_dir).parent.name,
                    file(params.flumut_reference_dir),
                ]
            )
            ch_target_reference = Channel.value(
                [
                    file(params.nextclade_datasets_dir).parent.name,
                    file(params.nextclade_datasets_dir),
                ]
            )
            coordinate_conversion(ch_source_reference.combine(ch_target_reference))

            // --- Step 8: Annotate each per-sample/per-segment VCF with FluMutDB markers ---
            //   Additive step: leaves the *_final.vcf produced above untouched and writes a
            //   separate *_flumut_merged.vcf carrying an extra FLUMUT= INFO field on rows whose
            //   AA change matches a curated FluMutDB marker. No rows are added or removed.
            //   Markers come from the flumutdb package, translated into Nextclade coordinates
            //   via the coordinate_conversion map, so the per-sample flumut TSV is not needed.
            //   Run per (sample, segment) so staging stays simple; each task consumes only its
            //   segment's slim coordinate-map shard.
            ch_vcf_flat = post_process_vcf_batch.out.vcf       // [segment, [vcfs]]
                .transpose()                                   // [segment, vcf]
                .map { segment, vcf ->
                    def sample_id = vcf.name.split('__')[0]
                    [sample_id, segment, vcf]
                }

            // The join below is an inner join on the segment label, and ch_shards only ever
            // carries the eight canonical labels. Split first so a VCF labelled anything else
            // (the 'unknown_segment' fallback of extract_segment) is reported rather than
            // vanishing, leaving a *_final.vcf with no merged counterpart and no explanation.
            ch_vcf_by_segment = ch_vcf_flat.branch { sample_id, segment, vcf ->
                known:   SEGMENTS.contains(segment)
                unknown: true
            }

            ch_vcf_by_segment.unknown.subscribe { sample_id, segment, vcf ->
                log.warn "merge_flumut_mutations: ${sample_id} has no coordinate map for segment " +
                         "'${segment}' (not one of ${SEGMENTS.join(', ')}); ${vcf.name} is left " +
                         "unannotated and produces no *_flumut_merged.vcf"
            }

            // Key each per-segment shard by its segment label parsed from the filename.
            ch_shards = coordinate_conversion.out.json_shards
                .flatten()
                .map { f -> [ (f.name =~ /coordinate_map\.([^.]+)\.json/)[0][1], f ] }  // [SEG, shard]

            ch_merge_in = ch_vcf_by_segment.known
                .map { sample_id, segment, vcf -> [segment, sample_id, vcf] }
                .combine(ch_shards, by: 0)                     // + coord shard (by segment)
                .map { segment, sample_id, vcf, shard -> [sample_id, segment, vcf, shard] }
                // -> [sample_id, segment, vcf, coord_shard]

            merge_flumut_mutations(ch_merge_in)

            // One merge task per (sample, segment), so group back to one entry per sample
            // before joining; it[2] is a list of YAMLs here, hence '+' rather than '<<'.
            ch_provenance = ch_provenance
                .join(merge_flumut_mutations.out.provenance.groupTuple())
                .map{ it -> [it[0], it[1] + it[2]] }

            // coordinate_conversion runs once per pipeline invocation, not per sample, so
            // its single provenance file is combined onto every sample rather than joined.
            ch_provenance = ch_provenance
                .combine(coordinate_conversion.out.provenance)
                .map{ it -> [it[0], it[1] << it[2]] }
        }

        collect_provenance(ch_provenance)

}
