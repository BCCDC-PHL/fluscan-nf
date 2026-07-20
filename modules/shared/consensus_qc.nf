// shared/consensus_qc.nf - Consensus FASTA quality control processes (shared by PAFTOOLS and NEXTCLADE workflows)

// prepare_consensus derives the authoritative isolate name from the FASTA
// header (the input sample_id is only a provisional, possibly-colliding, tag),
// so outputs are published under a directory named after that isolate, e.g.
//   out/WIN-AH-2022-FAV-1687-1-OS-2.fasta -> <outdir>/WIN-AH-2022-FAV-1687-1-OS-2/WIN-AH-2022-FAV-1687-1-OS-2.fasta
def publishUnderIsolate = { path ->
    def filename = file(path).name
    def isolate = filename.replaceAll(/(_qc_report)?\.(fasta|tsv)$/, '')
    "${isolate}/${filename}"
}

process prepare_consensus {

    tag { sample_id }

    errorStrategy 'ignore'

    conda "${projectDir}/environments/fluscan.yml"

    publishDir "${params.outdir}", pattern: "out/*.fasta",         mode: 'copy', saveAs: publishUnderIsolate
    publishDir "${params.outdir}", pattern: "out/*_qc_report.tsv", mode: 'copy', saveAs: publishUnderIsolate

    input:
    tuple val(sample_id), path(fasta_input)

    output:
    tuple env(isolate), path("out/*.fasta"),          emit: cleaned
    tuple env(isolate), path("out/*_qc_report.tsv"),  emit: qc
    tuple val(sample_id), path("${sample_id}_prepare_consensus_provenance.yml"), emit: provenance

    script:
    mask_degen = params.consensus_mask_degenerate ? '--mask-degenerate' : '--no-mask-degenerate'
    """
    printf -- "- process_name: prepare_consensus\\n"     > ${sample_id}_prepare_consensus_provenance.yml
    printf -- "  tools:\\n"                             >> ${sample_id}_prepare_consensus_provenance.yml
    printf -- "    - tool_name: prepare_consensus.py\\n" >> ${sample_id}_prepare_consensus_provenance.yml
    printf -- "      tool_version: 2.0.0\\n"            >> ${sample_id}_prepare_consensus_provenance.yml

    prepare_consensus.py \
        --input ${fasta_input} \
        --sample_id ${sample_id} \
        --outdir out \
        ${mask_degen} \
        --max-n-frac ${params.consensus_max_n_frac}

    isolate=\$(cat out/isolate_id.txt)
    """
}

process split_fasta {
    tag { sample_id }

    errorStrategy 'ignore'

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
