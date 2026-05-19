#!/usr/bin/env nextflow

/*
 *   Nextclade consensus FASTA workflow
 *   ------------------------------------
 *   Processes pre-assembled consensus FASTA files through quality control,
 *   per-segment Nextclade alignment against GooseGuangdong1996 reference
 *   datasets, and conversion of Nextclade JSON output to a freebayes-style VCF.
 *
 *   Expected input FASTA header format:
 *     >SAMPLE_ID_segmentN_SEGMENT
 *   e.g.
 *     >WIN-AH-2024-FAV-0415-1-OS_segment1_PB2
 *
 *   Segment datasets are resolved from --nextclade_datasets_dir by lowercasing
 *   the segment label (e.g. HA -> ha, M -> m).
 */

nextflow.enable.dsl = 2

include { prepare_consensus }                    from '../modules/consensus_qc.nf'
include { split_fasta as split_fasta_consensus } from '../modules/consensus_qc.nf'
include { flumut }                               from '../modules/flumut.nf'
include { nextclade_run }                        from '../modules/nextclade.nf'
include { nextclade_to_vcf }                     from '../modules/nextclade.nf'
include { post_process_vcf }                     from '../modules/variant_calling.nf'

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

        if (params.nextclade_datasets_dir == 'NO_FILE') {
            error "ERROR: --nextclade_datasets_dir is required for the Nextclade workflow."
        }

        // --- Step 1: Validate and clean input consensus FASTA ---
        prepare_consensus(ch_fasta_input)

        // --- Step 2: HPAI mutation screening (per sample, parallel branch) ---
        flumut(prepare_consensus.out.cleaned)

        // --- Step 3: Split cleaned FASTA into per-segment files ---
        ch_segment_with_dataset = split_fasta_consensus(prepare_consensus.out.cleaned)
            .transpose()
            .map(extract_segment)
            .map { sample_id, segment, fa ->
                def dataset_dir = file("${params.nextclade_datasets_dir}/${segment.toLowerCase()}")
                return [sample_id, segment, fa, dataset_dir]
            }

        // --- Step 4: Align each segment to its Nextclade dataset ---
        nextclade_run(ch_segment_with_dataset)

        // --- Step 5: Convert nextclade.json to freebayes-style VCF ---
        nextclade_to_vcf(nextclade_run.out.json_with_ref)

        // --- Step 6: Post-process VCF (e.g. fix MNV codons) ---
        post_process_vcf(
            nextclade_to_vcf.out.vcf.map{ sample_id, segment, vcf, dataset_dir ->
                def reference = file("${dataset_dir}/reference.fasta")
                return [sample_id, segment, vcf, reference]
            }
        )

}
