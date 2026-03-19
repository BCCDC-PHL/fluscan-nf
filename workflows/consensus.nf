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
 *     >accession|strain|subtype|segment
 *   e.g.
 *     >NC_026438.1|A/California/07/2009|H1N1|PB2
 *
 */

nextflow.enable.dsl = 2

include { pipeline_provenance } from '../modules/provenance.nf'
include { collect_provenance }  from '../modules/provenance.nf'
include { prepare_consensus }   from '../modules/consensus_qc.nf'
include { makeblastdb }         from '../modules/blast.nf'
include { blastn }              from '../modules/blast.nf'
include { consensus_paftools_call }   from '../modules/variant_calling.nf'
include { consensus_sort_vcf }        from '../modules/variant_calling.nf'
include { consensus_snpeff_annotation } from '../modules/variant_calling.nf'
include { consensus_reformat_vcf }    from '../modules/variant_calling.nf'


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

        if (params.snpeff_config == 'NO_FILE') {
            error "ERROR: --snpeff_config is required for the consensus workflow."
        }

        ch_snpeff_config = Channel.fromPath(params.snpeff_config, checkIfExists: true)

        // Provenance channel
        ch_provenance = ch_fasta_input.map{ it -> it[0] }

        // --- Step 1: Validate and clean input consensus FASTA ---
        prepare_consensus(ch_fasta_input)

        // --- Step 2: Build BLAST database (once, shared across all samples) ---
        makeblastdb(ch_nt_database)

        // --- Step 3: BLAST cleaned consensus against reference database ---
        blastn(prepare_consensus.out.cleaned, makeblastdb.out.db)

        // --- Step 4: Assign reference and snpEff DB from BLAST results ---
        ch_consensus_targets = prepare_consensus.out.cleaned
            .join(blastn.out.main.map { sample_id, blast_filter_tsv, reference_fasta ->
                [sample_id, blast_filter_tsv, reference_fasta]
            })
            .map { sample_id, cleaned_fasta, blast_filter_tsv, reference_fasta ->
                def rows = blast_filter_tsv.toFile().readLines()
                if (rows.size() < 2) {
                    throw new RuntimeException("ERROR [${sample_id}]: No filtered BLAST hits found in ${blast_filter_tsv}.")
                }

                def header = rows[0].split('\t')
                def sseqid_index = header.findIndexOf { it == 'sseqid' }
                def evalue_index = header.findIndexOf { it == 'evalue' }
                if (sseqid_index < 0 || evalue_index < 0) {
                    throw new RuntimeException("ERROR [${sample_id}]: BLAST filter TSV is missing required columns.")
                }

                def subtypes = rows.drop(1).collect { line ->
                    def fields = line.split('\t')
                    if (fields.size() <= Math.max(sseqid_index, evalue_index)) {
                        throw new RuntimeException("ERROR [${sample_id}]: Malformed BLAST filter row '${line}'.")
                    }

                    def evalue = fields[evalue_index] as Double
                    if (evalue > (params.consensus_blast_max_evalue as Double)) {
                        throw new RuntimeException("ERROR [${sample_id}]: BLAST hit exceeded max E-value (${fields[evalue_index]} > ${params.consensus_blast_max_evalue}).")
                    }

                    def sseqid_fields = fields[sseqid_index].split(/\|/)
                    if (sseqid_fields.size() < 4) {
                        throw new RuntimeException("ERROR [${sample_id}]: Unexpected BLAST subject header '${fields[sseqid_index]}'.")
                    }

                    sseqid_fields[-2].toLowerCase()
                }.unique()

                if (subtypes.size() != 1) {
                    throw new RuntimeException("ERROR [${sample_id}]: Expected exactly one subtype across BLAST hits, found ${subtypes}.")
                }

                def subtype = subtypes[0]
                def snpeff_db_name = "influenza_${subtype}"
                [sample_id, cleaned_fasta, reference_fasta, subtype, snpeff_db_name]
            }

        // --- Step 5: Align consensus to reference and call variants with paftools ---
        consensus_paftools_call(ch_consensus_targets)

        // --- Step 6: Sort VCF for downstream annotation ---
        consensus_sort_vcf(consensus_paftools_call.out.vcf)

        // --- Step 7: Annotate amino-acid impacts with SnpEff ---
        consensus_snpeff_annotation(consensus_sort_vcf.out.vcf.combine(ch_snpeff_config))

        // --- Step 8: Reformat annotated VCF to summary TSV ---
        consensus_reformat_vcf(consensus_snpeff_annotation.out.vcf)

        // --- Provenance collection ---
        ch_provenance = ch_provenance.combine(ch_pipeline_provenance).map{ it ->      [it[0], [it[1]]] }
        ch_provenance = ch_provenance.join(prepare_consensus.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
        ch_provenance = ch_provenance.join(blastn.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
        ch_provenance = ch_provenance.join(consensus_paftools_call.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
        ch_provenance = ch_provenance.join(consensus_sort_vcf.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
        ch_provenance = ch_provenance.join(consensus_snpeff_annotation.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
        ch_provenance = ch_provenance.join(consensus_reformat_vcf.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
        collect_provenance(ch_provenance)

}
