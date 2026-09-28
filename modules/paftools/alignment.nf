// paftools/alignment.nf - minimap2 alignment and PAF coverage processes for the PAFTOOLS workflow

process paftools_align {

    tag { sample_id }

    errorStrategy 'ignore'
    label 'heavy'

    input:
    tuple val(sample_id),
          path(consensus_fastas, stageAs: 'consensus/*'),
          path(reference_fastas, stageAs: 'ref_input/*')

    output:
    tuple val(sample_id), path('*.paf'), path('*fasta'), emit: paf
    tuple val(sample_id), path("${sample_id}_paftools_align_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: paftools_align\\n" > ${sample_id}_paftools_align_provenance.yml
    printf -- "  tools: \\n" >> ${sample_id}_paftools_align_provenance.yml
    printf -- "  - tool_name: minimap2\\n    tool_version: \$(minimap2 --version 2>&1 | head -n1)\\n" >> ${sample_id}_paftools_align_provenance.yml

    for consensus_fa in consensus/*.s.fasta; do
        consensus_basename=\$(basename "\$consensus_fa" .s.fasta)
        consensus_core="\${consensus_basename#${sample_id}_}"
        segment="\${consensus_core%_raw}"
        reference_fa="ref_input/${sample_id}_ref_\${segment}.s.fasta"
        output_prefix="${sample_id}_\${segment}"
        [ -s "\$reference_fa" ] || { echo "WARN: no reference for \$segment, skipping" >&2; continue; }
        minimap2 -cx asm20 --cs "\$reference_fa" "\$consensus_fa" > "\${output_prefix}_raw.paf" || { echo "WARN: minimap2 failed for \$segment" >&2; rm -f "\${output_prefix}_raw.paf"; continue; }
        sort -k6,6 -k8,8n "\${output_prefix}_raw.paf" > "\${output_prefix}.paf"
        rm "\${output_prefix}_raw.paf"
    done
    cp ref_input/*.fasta .
    """
}


process consensus_coverage_qc {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}/qc", mode: 'copy', pattern: "*.tsv"

    input:
    tuple val(sample_id),
          path(pafs, stageAs: 'pafs/*'),
          path(references, stageAs: 'ref_input/*')

    output:
    tuple val(sample_id), path("*_coverage_qc.tsv"), emit: coverage

    script:
    """
    for paf_file in pafs/*.paf; do
        paf_basename=\$(basename "\$paf_file" .paf)
        segment="\${paf_basename#${sample_id}_}"
        reference_fa="ref_input/${sample_id}_ref_\${segment}.s.fasta"
        consensus_coverage_qc.py \\
          --paf "\$paf_file" \\
          --reference "\$reference_fa" \\
          --sample-id ${sample_id} \\
          --segment "\$segment" \\
          --pass-pct ${params.consensus_coverage_pass_pct} \\
          --warn-pct ${params.consensus_coverage_warn_pct} \\
          --out-tsv ${sample_id}_\${segment}_coverage_qc.tsv
    done
    """
}
