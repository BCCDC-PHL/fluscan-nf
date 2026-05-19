// shared/consensus_qc.nf - Consensus FASTA quality control processes (shared by PAFTOOLS and NEXTCLADE workflows)

process prepare_consensus {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}/", pattern: "${sample_id}_cleaned.fasta", mode: 'copy'
    publishDir "${params.outdir}/${sample_id}/", pattern: "${sample_id}_qc_report.tsv", mode: 'copy'

    input:
    tuple val(sample_id), path(fasta_input)

    output:
    tuple val(sample_id), path("${sample_id}_cleaned.fasta"), emit: cleaned
    tuple val(sample_id), path("${sample_id}_qc_report.tsv"), emit: qc
    tuple val(sample_id), path("${sample_id}_prepare_consensus_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: prepare_consensus\\n"     > ${sample_id}_prepare_consensus_provenance.yml
    printf -- "  tools:\\n"                             >> ${sample_id}_prepare_consensus_provenance.yml
    printf -- "    - tool_name: prepare_consensus.py\\n" >> ${sample_id}_prepare_consensus_provenance.yml
    printf -- "      tool_version: 1.0.0\\n"            >> ${sample_id}_prepare_consensus_provenance.yml

    mv ${fasta_input} raw.fasta

    prepare_consensus.py \
        --input raw.fasta \
        --sample_id ${sample_id} \
        --output ${sample_id}_cleaned.fasta \
        --output_qc ${sample_id}_qc_report.tsv
    """
}

process split_fasta {
    tag { sample_id }

    //publishDir "${params.outdir}/${sample_id}/", pattern: "${sample_id}_s_*.fasta", mode: 'copy'

    input:
    tuple val(sample_id), path(fasta)

    output:
    tuple val(sample_id), path("${sample_id}*.s.fasta"), emit: fastas

    script:
    """
    split_fasta.py split ${fasta}
    """
}
