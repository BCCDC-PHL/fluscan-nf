// flumut.nf - FluMut influenza molecular marker detection

process flumut {

    tag { sample_id }

    conda "${projectDir}/environments/fluscan.yml"

    publishDir "${params.outdir}/${sample_id}/flumut", pattern: "${sample_id}_flumut_*.tsv", mode: 'copy'

    input:
    tuple val(sample_id), path(consensus_fasta)

    output:
    tuple val(sample_id), path("${sample_id}_flumut_markers.tsv"),   emit: markers
    tuple val(sample_id), path("${sample_id}_flumut_mutations.tsv"),  emit: mutations
    tuple val(sample_id), path("${sample_id}_flumut_wide.tsv"),       emit: wide_mutations
    tuple val(sample_id), path("${sample_id}_flumut_literature.tsv"), emit: literature
    tuple val(sample_id), path("${sample_id}_flumut_provenance.yml"), emit: provenance

    script:
    def regex_opt = task.ext.segment_regex ? "--name-regex '${task.ext.segment_regex}'" : ""
    """
    printf -- "- process_name: flumut\\n"                                                        > ${sample_id}_flumut_provenance.yml
    printf -- "  tools:\\n"                                                                      >> ${sample_id}_flumut_provenance.yml
    printf -- "    - tool_name: flumut\\n"                                                       >> ${sample_id}_flumut_provenance.yml
    printf -- "      tool_version: \$(flumut --version 2>&1 | cut -d' ' -f3)\\n"                >> ${sample_id}_flumut_provenance.yml

    sed -E -i 's/([\\|_])M\$/\\1MP/g' ${consensus_fasta}

    flumut \\
        -m ${sample_id}_flumut_markers.tsv \\
        -M ${sample_id}_flumut_wide.tsv \\
        -l ${sample_id}_flumut_literature.tsv \\
        --skip-unmatch-names \\
        --skip-unknown-segments ${regex_opt} \\
        ${consensus_fasta}

    melt_flumut.py --input ${sample_id}_flumut_wide.tsv --output ${sample_id}_flumut_mutations.tsv
    """
}
