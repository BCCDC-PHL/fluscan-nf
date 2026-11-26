// shared/genoflu.nf - GenoFLU H5N1 genotype assignment (shared by all workflows)

process genoflu {

    tag { sample_id }

    conda "${projectDir}/environments/genoflu.yml"

    publishDir "${params.outdir}/${sample_id}/genoflu", pattern: "${sample_id}_genoflu.tsv", mode: 'copy'

    input:
    tuple val(sample_id), path(consensus_fasta)

    output:
    tuple val(sample_id), path("${sample_id}_genoflu.tsv"),         emit: tsv
    tuple val(sample_id), path("${sample_id}_genoflu_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: genoflu\\n"                                                       > ${sample_id}_genoflu_provenance.yml
    printf -- "  tools:\\n"                                                                      >> ${sample_id}_genoflu_provenance.yml
    printf -- "    - tool_name: genoflu\\n"                                                      >> ${sample_id}_genoflu_provenance.yml
    printf -- "      tool_version: \$(genoflu.py --version 2>&1 | awk '{print \$NF}')\\n"        >> ${sample_id}_genoflu_provenance.yml

    genoflu.py -f ${consensus_fasta} -n ${sample_id}

    mv ${sample_id}_*_stats.tsv ${sample_id}_genoflu.tsv
    """
}
