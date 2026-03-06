#!/usr/bin/env nextflow

/*
 *   Consensus FASTA analysis workflow
 *   ------------------------------------
 *   Processes pre-assembled consensus FASTA files through quality control,
 *   reference selection via BLAST, and clade assignment via Nextclade.
 *
 *   Expected input FASTA header format:
 *     >SAMPLE_ID_segmentN_SEGMENT
 *   e.g.
 *     >WIN-AH-2024-FAV-0415-1-OS_segment1_PB2
 *     >R3980461859-3370-A-D05_segment1_PB2
 *
 *   Reference database header format (--nt_database):
 *     >accession|strain|subtype|segment
 *   e.g.
 *     >NC_026438.1|A/California/07/2009|H1N1|PB2
 *
 *   Nextclade datasets directory (--nextclade_datasets_dir):
 *     <datasets_dir>/<subtype_lower>/<segment_lower>/
 *   e.g.
 *     h1n1/ha  h1n1/mp  h1n1/pb2  ...
 *     h3n2/ha  h3n2/mp  h3n2/pb2  ...
 */

nextflow.enable.dsl = 2

include { pipeline_provenance } from '../modules/provenance.nf'
include { collect_provenance }  from '../modules/provenance.nf'
include { prepare_consensus }   from '../modules/consensus_qc.nf'
include { nextclade_consensus } from '../modules/consensus_qc.nf'
include { makeblastdb }         from '../modules/blast.nf'
include { blastn }              from '../modules/blast.nf'


workflow CONSENSUS {

    ch_workflow_metadata = Channel.value([
        workflow.sessionId,
        workflow.runName,
        workflow.manifest.name,
        workflow.manifest.version,
        workflow.start,
    ])

    ch_pipeline_provenance = pipeline_provenance(ch_workflow_metadata)

    // Discover consensus FASTA files in the input directory.
    // Each file is expected to contain 8 segments for a single sample.
    ch_fasta_input = Channel.fromPath(
        params.fasta_input + "/*.{fasta,fa,fna}",
        checkIfExists: true
    ).map { fasta ->
        def sample_id = fasta.baseName.replaceAll(/_consensus$/, '')
        return [sample_id, fasta]
    }

    ch_nt_database = Channel.fromPath(params.nt_database)

    main:

        // Provenance channel
        ch_provenance = ch_fasta_input.map{ it -> it[0] }

        // --- Step 1: Validate and clean input consensus FASTA ---
        prepare_consensus(ch_fasta_input)

        // --- Step 2: Build BLAST database (once, shared across all samples) ---
        makeblastdb(ch_nt_database)

        // --- Step 3: BLAST cleaned consensus against reference database ---
        //     The filter TSV will contain the best reference for each segment,
        //     including the subtype and segment label in the reference header.
        blastn(prepare_consensus.out.cleaned, makeblastdb.out.db)

        // --- Step 4: Run Nextclade for each segment ---
        //     Combines cleaned FASTA with the filtered BLAST TSV to determine
        //     the appropriate Nextclade dataset (subtype/segment) per sequence.
        //     blastn.out.main emits: [sample_id, filter_tsv, ref_fasta]
        ch_blast_filter_tsv = blastn.out.main.map { it -> [it[0], it[1]] }
        ch_for_nextclade = prepare_consensus.out.cleaned.join(ch_blast_filter_tsv)

        nextclade_consensus(
            ch_for_nextclade,
            params.nextclade_datasets_dir
        )

        // --- Provenance collection ---
        ch_provenance = ch_provenance.combine(ch_pipeline_provenance).map{ it ->      [it[0], [it[1]]] }
        ch_provenance = ch_provenance.join(prepare_consensus.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
        ch_provenance = ch_provenance.join(nextclade_consensus.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
        collect_provenance(ch_provenance)

}
