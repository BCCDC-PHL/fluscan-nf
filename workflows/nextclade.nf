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
include { consensus_identity_qc }                from '../modules/shared/consensus_qc.nf'
include { flumut }                               from '../modules/shared/flumut.nf'
include { nextclade_run_batch }                  from '../modules/nextclade/nextclade.nf'
include { nextclade_to_vcf_batch }               from '../modules/nextclade/nextclade.nf'
include { post_process_vcf_batch }               from '../modules/shared/post_process.nf'
include { coordinate_conversion }                from '../modules/nextclade/coordinate_conversion.nf'
include { insert_flumut_mutations }              from '../modules/nextclade/insert_flumut_mutations.nf'

workflow NEXTCLADE {

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

        if (
            !params.nextclade_datasets_dir ||
            params.nextclade_datasets_dir == 'NO_FILE' ||
            !file(params.nextclade_datasets_dir).exists()
        ) {
            error "ERROR: --nextclade_datasets_dir must be set to an existing path for the Nextclade workflow."
        }

        if (
            !params.flumut_reference_dir ||
            params.flumut_reference_dir == 'NO_FILE' ||
            !file(params.flumut_reference_dir).exists()
        ) {
            error "ERROR: --flumut_reference_dir must be set to an existing path (flumut's Nextclade dataset directory, the source coordinate reference) for the Nextclade workflow."
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

        // --- Step 8: Insert flumut AA markers into each per-sample/per-segment VCF ---
        //   New, additive step: leaves the *_final.vcf produced above untouched and writes a
        //   separate *_flumut_merged.vcf. Labels every AA change DETECTED_BY=fluscan|flumut|both
        //   and appends flumut-only markers, remapped into Nextclade coordinates via the
        //   coordinate_conversion map. Run per (sample, segment) so staging stays simple.
        //   Each task consumes only its segment's slim coordinate-map shard (target CDS
        //   placement is carried in the shard, so no GFF/dataset dir is needed here).
        ch_vcf_flat = post_process_vcf_batch.out.vcf       // [segment, [vcfs], dataset_dir]
            .transpose()                                   // [segment, vcf, dataset_dir]
            .map { segment, vcf, dataset_dir ->
                def sample_id = vcf.name.split('__')[0]
                [sample_id, segment, vcf]                  // dataset_dir no longer needed
            }

        // Key each per-segment shard by its segment label parsed from the filename.
        ch_shards = coordinate_conversion.out.json_shards
            .flatten()
            .map { f -> [ (f.name =~ /coordinate_map\.([^.]+)\.json/)[0][1], f ] }  // [SEG, shard]

        ch_insert_in = ch_vcf_flat
            .map { sample_id, segment, vcf -> [segment, sample_id, vcf] }
            .combine(ch_shards, by: 0)                     // + coord shard (by segment)
            .map { segment, sample_id, vcf, shard -> [sample_id, segment, vcf, shard] }
            .combine(flumut.out.mutations, by: 0)          // + flumut TSV (by sample)
            .combine(prepare_consensus.out.cleaned, by: 0) // + consensus FASTA (by sample)
            // -> [sample_id, segment, vcf, coord_shard, mutations, consensus]

        insert_flumut_mutations(ch_insert_in)

}
