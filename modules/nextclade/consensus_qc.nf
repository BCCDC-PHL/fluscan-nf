// nextclade/consensus_qc.nf - Consensus/reference identity QC for the NEXTCLADE workflow

// Report-only consensus/reference identity QC for a single segment.
process consensus_identity_qc {
    tag { "${sample_id}|${segment}" }

    errorStrategy 'ignore'

    conda "${projectDir}/environments/fluscan.yml"

    publishDir "${params.outdir}/${sample_id}/qc", mode: 'copy', pattern: "*.tsv"

    input:
    tuple val(sample_id), val(segment), path(consensus_fasta), path(reference_fasta)

    output:
    tuple val(sample_id), val(segment), path("*_identity_qc.tsv"), emit: tsv

    script:
    """
    consensus_identity_qc.py \
        --consensus ${consensus_fasta} \
        --reference ${reference_fasta} \
        --sample-id ${sample_id} \
        --segment ${segment} \
        --pass-pct ${params.consensus_identity_pass_pct} \
        --warn-pct ${params.consensus_identity_warn_pct} \
        --out-tsv ${sample_id}_${segment}_identity_qc.tsv
    """
}
