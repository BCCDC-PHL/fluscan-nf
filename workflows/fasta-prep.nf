#!/usr/bin/env nextflow

/*
 *   Consensus FASTA preparation, shared by PAFTOOLS and NEXTCLADE
 *   -------------------------------------------------------------
 *   Hashes each input FASTA, validates and cleans it, aggregates the consensus QC rows and
 *   runs FluMut, once per sample however many FASTA branches consume the result.
 *
 *   The sample name comes from the input filename (modules/shared/sample_id.nf); input
 *   headers supply only the segment label. Accepted header formats:
 *     Underscore-delimited (legacy): >SAMPLE_ID_segmentN_SEGMENT
 *       e.g. >WIN-AH-2024-FAV-0415-1-OS_segment1_PB2
 *     Pipe-delimited (GISAID-style): >EPI_ISL_ACCESSION|subtype|SEGMENT|...|clade_SEGMENT
 *       e.g. >EPI_ISL_20399898|A_/_H5N1|NA|6|EPI5261179|2.3.4.4b_NA
 *   Cleaned headers are normalised to {sample_id}_{SEGMENT} by prepare_consensus.
 */

include { hash_files }              from '../modules/shared/hash_files.nf'
include { prepare_consensus }       from '../modules/shared/consensus_qc.nf'
include { flumut }                  from '../modules/shared/flumut.nf'
include { consensus_fasta_channel } from '../modules/shared/sample_id.nf'
include { samplesheet_channel }     from '../modules/shared/sample_id.nf'

workflow FASTA_PREP {

    main:

        // Consensus FASTAs from the sample sheet or the input directory.
        // Each file is expected to contain 8 segments for a single sample.
        ch_fasta_input = params.samplesheet_input != 'NO_FILE'
            ? samplesheet_channel(params.samplesheet_input, ['FASTA'])
            : consensus_fasta_channel(params.fasta_input)

        // Hash the input FASTA. One file per sample, so no list-wrapping here.
        hash_files(ch_fasta_input.combine(Channel.of("fasta_input")))

        // Validate and clean input consensus FASTA.
        prepare_consensus(ch_fasta_input)

        // One run-level table of the per-segment consensus QC rows.
        prepare_consensus.out.qc
            .map{ it -> it[1] }
            .collectFile(
                name: 'consensus_qc_summary.tsv',
                storeDir: "${params.outdir}/aggregate",
                keepHeader: true,                   // keep header from first file only
                skip: 1,                            // drop the header line from the rest
                sort: false
            )

        // HPAI mutation screening (per sample).
        flumut(prepare_consensus.out.cleaned)

        // Inner joins: prepare_consensus runs errorStrategy 'ignore', so a sample it rejects
        // leaves every branch's provenance chain.
        ch_provenance = hash_files.out.provenance
            .join(prepare_consensus.out.provenance)
            .join(flumut.out.provenance)
            .map { sample_id, hash, prepare, mutations -> [sample_id, [hash, prepare, mutations]] }

    emit:
        cleaned    = prepare_consensus.out.cleaned    // [sample_id, cleaned consensus FASTA]
        provenance = ch_provenance                    // [sample_id, [hash, prepare, flumut]]
}
