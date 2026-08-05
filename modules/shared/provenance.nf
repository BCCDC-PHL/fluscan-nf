// shared/provenance.nf - Pipeline and sample provenance tracking (shared by all workflows)

process collect_provenance {

  tag { sample_id }

  executor 'local'

  publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_*_provenance.yml", mode: 'copy'

  input:
  tuple val(sample_id), path(provenance_files)

  output:
  tuple val(sample_id), file("${sample_id}_*_provenance.yml")

  script:
  """
  cat ${provenance_files} > ${sample_id}_\$(date +%Y%m%d%H%M%S)_provenance.yml
  """
}

process pipeline_provenance {

    tag { pipeline_name + " / " + pipeline_version }

    executor 'local'

    input:
    tuple val(session_id), val(run_name), val(pipeline_name), val(pipeline_version), val(timestamp_analysis_start), val(databases)

    output:
    file("pipeline_provenance.yml")

    script:
    // `databases` is an insertion-ordered map of reference-database facts, so which keys
    // are recorded follows the workflow: the flumut bundle only appears when the merge ran.
    // Values are single-quoted so a date-like version ("2025-07-07") reads back as the
    // string it is rather than as a YAML date.
    def database_lines = databases.collect { key, value ->
        "    printf -- \"  ${key}: '${value}'\\\\n\" >> pipeline_provenance.yml"
    }.join('\n')
    """
    printf -- "- pipeline_name: ${pipeline_name}\\n"                       >> pipeline_provenance.yml
    printf -- "  pipeline_version: ${pipeline_version}\\n"                 >> pipeline_provenance.yml
    printf -- "  nextflow_session_id: ${session_id}\\n"                    >> pipeline_provenance.yml
    printf -- "  nextflow_run_name: ${run_name}\\n"                        >> pipeline_provenance.yml
    printf -- "  timestamp_analysis_start: ${timestamp_analysis_start}\\n" >> pipeline_provenance.yml
${database_lines}
    """
}
