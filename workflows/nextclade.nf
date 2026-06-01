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
include { flumut }                               from '../modules/shared/flumut.nf'
include { nextclade_run_batch }                  from '../modules/nextclade/nextclade.nf'
include { nextclade_to_vcf_batch }               from '../modules/nextclade/nextclade.nf'
include { post_process_vcf_batch }               from '../modules/shared/post_process.nf'

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

        // --- Step 1: Validate and clean input consensus FASTA ---
        prepare_consensus(ch_fasta_input)

        // --- Step 2: HPAI mutation screening (per sample, parallel branch) ---
        flumut(prepare_consensus.out.cleaned)

        // --- Step 3: Split cleaned FASTA into per-segment files, then group by
        //             segment across all samples for cross-sample batching ---
        ch_segment_batches = split_fasta_consensus(prepare_consensus.out.cleaned)
            .transpose()
            .map(extract_segment)                          // [sample_id, segment, fa]
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

}
