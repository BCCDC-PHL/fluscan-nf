process cutadapt {

    tag { sample_id }

    input:
    tuple val(sample_id), path(reads_1), path(reads_2), path(primers)

    output:
    tuple val(sample_id), path("${sample_id}_R1.trim_p.fastq.gz"), path("${sample_id}_R2.trim_p.fastq.gz"), emit: primer_trimmed_reads
    path("${sample_id}.cutadapt.log"), emit: log
    tuple val(sample_id), path("${sample_id}_cutadapt_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: cutadapt\\n"    >> ${sample_id}_cutadapt_provenance.yml
    printf -- "  tools:\\n"                    >> ${sample_id}_cutadapt_provenance.yml
    printf -- "    - tool_name: cutadapt\\n"   >> ${sample_id}_cutadapt_provenance.yml
    printf -- "      tool_version: \$(cutadapt --version 2>&1 | cut -d ' ' -f 2)\\n" >> ${sample_id}_cutadapt_provenance.yml

    cutadapt \
      -j ${task.cpus} \
      -a file:${params.primers} \
      -A file:${params.primers} \
      -g file:${params.primers} \
      -G file:${params.primers} \
      -o ${sample_id}_R1.trim_p.fastq.gz \
      -p ${sample_id}_R2.trim_p.fastq.gz \
      ${sample_id}_R1.trim.fastq.gz \
      ${sample_id}_R2.trim.fastq.gz \
      > ${sample_id}.cutadapt.log
    """
}

process fastp {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "*fastp.*",  mode:'copy'

    input:
    tuple val(sample_id), path(reads_1), path(reads_2)

    output:
    tuple val(sample_id), path("${sample_id}_R1.trim.fastq.gz"), path("${sample_id}_R2.trim.fastq.gz"), emit: trimmed_reads
    path("${sample_id}.fastp.json"), emit: json
    path("${sample_id}.fastp.html"), emit: html
    tuple val(sample_id), path("${sample_id}_fastp_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: fastp\\n"  >> ${sample_id}_fastp_provenance.yml
    printf -- "  tools:\\n"               >> ${sample_id}_fastp_provenance.yml
    printf -- "    - tool_name: fastp\\n" >> ${sample_id}_fastp_provenance.yml
    printf -- "      tool_version: \$(fastp --version 2>&1 | cut -d ' ' -f 2)\\n" >> ${sample_id}_fastp_provenance.yml

    fastp \
      -t ${task.cpus} \
      -i ${reads_1} \
      -I ${reads_2} \
      -o ${sample_id}_R1.trim.fastq.gz \
      -O ${sample_id}_R2.trim.fastq.gz \
      -j ${sample_id}.fastp.json \
      -h ${sample_id}.fastp.html \
      --detect_adapter_for_pe \
      --cut_tail \
      --report_title "${sample_id} fastp report" 
    """
}


process fastp_json_to_csv {

  tag { sample_id }

  executor 'local'

  publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_fastp.csv", mode: 'copy'

  input:
  tuple val(sample_id), path(fastp_json)

  output:
  tuple val(sample_id), path("${sample_id}_fastp.csv")

  script:
  """
  fastp_json_to_csv.py -s ${sample_id} ${fastp_json} > ${sample_id}_fastp.csv
  """
}


process fastqc {

    tag { sample_id }

    input:
    tuple val(sample_id), path(reads_1), path(reads_2)

    output:
    //tuple val(sample_id), path("${sample_id}_R1.trim_p_fastqc.html"), path("${sample_id}_R2.trim_p_fastqc.html"), emit: html
    //tuple val(sample_id), path("${sample_id}_R1.trim_p_fastqc.zip"), path("${sample_id}_R2.trim_p_fastqc.zip"), emit: zip
    tuple val(sample_id), path("${sample_id}_R1.trim_p_fastqc.zip"), path("${sample_id}_R2.trim_p_fastqc.zip"), emit: zip
    tuple val(sample_id), path("${sample_id}_fastqc_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: fastqc\\n"  >> ${sample_id}_fastqc_provenance.yml
    printf -- "  tools:\\n"                >> ${sample_id}_fastqc_provenance.yml
    printf -- "    - tool_name: fastqc\\n" >> ${sample_id}_fastqc_provenance.yml
    printf -- "      tool_version: \$(fastqc --version 2>&1 | sed -n '1 p')\\n" >> ${sample_id}_fastqc_provenance.yml

    mkdir -p ./tmp

    fastqc \
      --threads ${task.cpus} \
      --dir ./tmp \
      ${sample_id}_R1.trim_p.fastq.gz \
      ${sample_id}_R2.trim_p.fastq.gz 
    
    """
}


process multiqc {

    publishDir "${params.outdir}", pattern: "*_multiqc_report.html", mode:'copy'

    input:
    path '*'

    output:
    path '*_multiqc_*'

    script:
    """
    multiqc . -n ${params.run_name}_fluviewer-nf_multiqc_report.html
    """
}
