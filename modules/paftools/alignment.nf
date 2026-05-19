// paftools/alignment.nf - minimap2 alignment and PAF coverage processes for the PAFTOOLS workflow

process paftools_align {

    tag { "${name}" }

    errorStrategy 'ignore'
    label 'heavy'

    input:
    tuple val(sample_id), val(segment), path(consensus_fasta), path(reference_fasta)

    output:
    tuple val(sample_id), val(segment), path("${name}.paf"), path(reference_fasta), emit: paf
    tuple val(sample_id), path("${provenance_file}"), emit: provenance

    script:
    provenance_file = "${sample_id}_paftools_align_provenance.yml"
    name = "${sample_id}_${segment}"
    """
    printf -- "- process_name: paftools_align\\n" > ${provenance_file}
    printf -- "  tools: \\n" >> ${provenance_file}
    printf -- "  - tool_name: minimap2\\n    tool_version: \$(minimap2 --version 2>&1 | head -n1)\\n" >> ${provenance_file}

    minimap2 -cx asm20 --cs ${reference_fasta} ${consensus_fasta} > ${name}_raw.paf
    sort -k6,6 -k8,8n ${name}_raw.paf > ${name}.paf
    rm ${name}_raw.paf

    """
}


process paftools_coverage_check {

    tag { "${sample_id}__${segment}" }

    publishDir "${params.outdir}/${sample_id}/alignment_qc", mode: 'copy', pattern: "*.tsv"

    input:
    tuple val(sample_id), val(segment), path(paf), path(reference)

    output:
    tuple val(sample_id), val(segment), path("${sample_id}_${segment}_coverage.tsv"), emit: coverage

    script:
    """
    paf_coverage_check.py \\
      --paf ${paf} \\
      --sample ${sample_id} \\
      --segment ${segment} \\
      --reference ${reference} \\
      --pass-thresh ${params.min_segment_coverage} \\
      --warn-thresh ${params.min_segment_coverage_strict} \\
      --out-tsv ${sample_id}_${segment}_coverage.tsv
    """
}


process paftools_aggregate_coverage_report {

    publishDir "${params.outdir}/_run_summary", mode: 'copy'

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
