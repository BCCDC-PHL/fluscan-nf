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
        segment=\$(basename "\$consensus_fa" .s.fasta | sed -E 's/.*_([^_]+)\$/\\1/')
        reference_fa="ref_input/${sample_id}_\${segment}.s.fasta"
        output_prefix="${sample_id}_\${segment}"
        minimap2 -cx asm20 --cs "\$reference_fa" "\$consensus_fa" > "\${output_prefix}_raw.paf" || { echo "WARN: minimap2 failed for \$segment" >&2; continue; }
        sort -k6,6 -k8,8n "\${output_prefix}_raw.paf" > "\${output_prefix}.paf"
        rm "\${output_prefix}_raw.paf"
    done
    cp ref_input/*.fasta .
    """
}


process paftools_coverage_check {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}/alignment_qc", mode: 'copy', pattern: "*.tsv"

    input:
    tuple val(sample_id),
          path(pafs, stageAs: 'pafs/*'),
          path(references, stageAs: 'ref_input/*')

    output:
    tuple val(sample_id), path("*_coverage.tsv"), emit: coverage

    script:
    """
    for paf_file in pafs/*.paf; do
        segment=\$(basename "\$paf_file" .paf | sed -E 's/.*_([^_]+)\$/\\1/')
        reference_fa="ref_input/${sample_id}_\${segment}.s.fasta"
        paf_coverage_check.py \\
          --paf "\$paf_file" \\
          --sample ${sample_id} \\
          --segment "\$segment" \\
          --reference "\$reference_fa" \\
          --pass-thresh ${params.min_segment_coverage} \\
          --warn-thresh ${params.min_segment_coverage_strict} \\
          --out-tsv ${sample_id}_\${segment}_coverage.tsv
    done
    """
}


process paftools_aggregate_coverage_report {

    publishDir "${params.outdir}/aggregate", mode: 'copy'

    input:
    path(tsv_files)

    output:
    path("run_coverage_summary.tsv"), emit: tsv
    path("run_coverage_summary.md"),  emit: md

    script:
    """
    aggregate_coverage_report.py \\
      --out-tsv run_coverage_summary.tsv \\
      --out-md  run_coverage_summary.md \\
      ${tsv_files}
    """
}
