// shared/consensus_qc.nf - Consensus FASTA quality control processes (shared by PAFTOOLS and NEXTCLADE workflows)

// sample_id arrives final from modules/shared/sample_id.nf and names everything
// this process emits, so the cleaned FASTA carries {sample_id}_{SEGMENT} headers
// and publishes to <outdir>/<sample_id>/. The input is staged under input/ so a
// file already called ${sample_id}.fasta cannot collide with the output.
process prepare_consensus {

    tag { sample_id }

    errorStrategy 'ignore'

    conda "${projectDir}/environments/fluscan.yml"

    publishDir "${params.outdir}/${sample_id}", pattern: "*.fasta",         mode: 'copy'
    publishDir "${params.outdir}/${sample_id}/qc", pattern: "*_qc_report.tsv", mode: 'copy'

    input:
    tuple val(sample_id), path(fasta_input, stageAs: 'input/*')

    output:
    tuple val(sample_id), path("${sample_id}.fasta"),         emit: cleaned
    tuple val(sample_id), path("${sample_id}_qc_report.tsv"), emit: qc
    tuple val(sample_id), path("${sample_id}_prepare_consensus_provenance.yml"), emit: provenance

    script:
    mask_degen = params.consensus_mask_degenerate ? '--mask-degenerate' : '--no-mask-degenerate'
    """
    printf -- "- process_name: prepare_consensus\\n"     > ${sample_id}_prepare_consensus_provenance.yml
    printf -- "  tools:\\n"                             >> ${sample_id}_prepare_consensus_provenance.yml
    printf -- "    - tool_name: prepare_consensus.py\\n" >> ${sample_id}_prepare_consensus_provenance.yml
    printf -- "      tool_version: 3.1.0\\n"            >> ${sample_id}_prepare_consensus_provenance.yml

    prepare_consensus.py \
        --input ${fasta_input} \
        --sample-id ${sample_id} \
        --output-fasta ${sample_id}.fasta \
        --output-qc ${sample_id}_qc_report.tsv \
        ${mask_degen} \
        --max-n-frac ${params.consensus_max_n_frac}
    """
}

process split_fasta {
    tag { sample_id }

    errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(fasta)

    output:
    tuple val(sample_id), path("${sample_id}*.s.fasta"), emit: fastas

    script:
    """
    split_fasta.py split ${fasta}
    """
}
