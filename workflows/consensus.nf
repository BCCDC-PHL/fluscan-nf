#!/usr/bin/env nextflow

/*
 *   Consensus FASTA analysis workflow
 *   ------------------------------------
 *   Processes pre-assembled consensus FASTA files through quality control,
 *   subtype-aware reference selection via BLAST, minimap2/paftools variant
 *   calling, bcftools VCF sorting, SnpEff annotation, and summary extraction.
 *
 *   Expected input FASTA header format:
 *     >SAMPLE_ID_segmentN_SEGMENT
 *   e.g.
 *     >WIN-AH-2024-FAV-0415-1-OS_segment1_PB2
 *     >R3980461859-3370-A-D05_segment1_PB2
 *
 *   Reference database header format (--nt_database):
 *     >accession|subtype|segment
 *   e.g.
 *     >NC_026438.1|H1N1|PB2
 *
 */

nextflow.enable.dsl = 2

include { pipeline_provenance } from '../modules/provenance.nf'
include { collect_provenance }  from '../modules/provenance.nf'
include { prepare_consensus }   from '../modules/consensus_qc.nf'
include { split_fasta as split_fasta_reference }         from '../modules/consensus_qc.nf'
include { split_fasta as split_fasta_consensus }         from '../modules/consensus_qc.nf'
include { map_reads_long }      from '../modules/alignment.nf'
include { makeblastdb }         from '../modules/blast.nf'
include { blastn }              from '../modules/blast.nf'
include { genoflu }              from '../modules/genoflu.nf'
include { flumut }               from '../modules/flumut.nf'
include { consensus_paftools_call }   from '../modules/variant_calling.nf'
include { consensus_sort_vcf }        from '../modules/variant_calling.nf'
include { consensus_snpeff_annotation } from '../modules/variant_calling.nf'
include { reformat_vcf }               from '../modules/variant_calling.nf'


workflow CONSENSUS {

    def extractSegment = { sample_id, fasta ->
        def segment = fasta =~ /.*_([A-z0-9-]+)\.s.*/
        segment = segment[0][1] ? segment[0][1] : "unknown_segment"
        return [sample_id, segment, fasta]
    }

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
        def sample_id = fasta.baseName.replaceAll(/[_\.]consensus$/, '').replaceAll(/_/,'-')
        sample_id = sample_id.replaceFirst(/^\d{2}-\d+(?:-IP\d+)?-/, "")
        return [sample_id, fasta]
    }

    ch_nt_database = Channel.fromPath(params.nt_database)

    main:

        if (params.snpeff_config == 'NO_FILE') {
            error "ERROR: --snpeff_config is required for the consensus workflow."
        }

        ch_snpeff_config = Channel.fromPath(params.snpeff_config, checkIfExists: true)

        // Provenance channel
        ch_provenance = ch_fasta_input.map{ it -> it[0] }

        // --- Step 1: Validate and clean input consensus FASTA ---
        prepare_consensus(ch_fasta_input)

        // --- Step 2: Build BLAST database (once, shared across all samples) ---
        makeblastdb(ch_nt_database).first().set{ch_blastdb}

        // --- Step 3: BLAST cleaned consensus against reference database ---
        blastn(prepare_consensus.out.cleaned, ch_blastdb)

        genoflu(prepare_consensus.out.cleaned)

        flumut(prepare_consensus.out.cleaned)

        ch_consensus_segment = split_fasta_consensus(prepare_consensus.out.cleaned).transpose().map(extractSegment)
        ch_reference_segment = split_fasta_reference(blastn.out.ref).transpose().map(extractSegment)
        
        ch_segment_sequences = ch_consensus_segment.join(ch_reference_segment, by: [0,1]) 

        map_reads_long(ch_segment_sequences)

        // --- Step 5: Align consensus to reference and call variants with paftools ---
        consensus_paftools_call(map_reads_long.out.paf)

        // --- Step 6: Sort VCF for downstream annotation ---
        consensus_sort_vcf(consensus_paftools_call.out.vcf)

        // --- Step 7: Annotate amino-acid impacts with SnpEff ---
        consensus_snpeff_annotation(consensus_sort_vcf.out.vcf.combine(ch_snpeff_config))

        // // --- Step 8: Reformat annotated VCF to summary TSV ---
        //reformat_vcf(consensus_snpeff_annotation.out.vcf)

        // // --- Provenance collection ---
        // ch_provenance = ch_provenance.combine(ch_pipeline_provenance).map{ it ->      [it[0], [it[1]]] }
        // ch_provenance = ch_provenance.join(prepare_consensus.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
        // ch_provenance = ch_provenance.join(blastn.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
        // ch_provenance = ch_provenance.join(consensus_paftools_call.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
        // ch_provenance = ch_provenance.join(consensus_sort_vcf.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
        // ch_provenance = ch_provenance.join(consensus_snpeff_annotation.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
        // ch_provenance = ch_provenance.join(reformat_vcf.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
        // collect_provenance(ch_provenance)

}
