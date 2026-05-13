// nextclade.nf - Per-segment Nextclade alignment and VCF conversion processes

process nextclade_run {

    conda "${projectDir}/environments/nextclade.yml"

    tag { "${sample_id}_${segment}" }

    publishDir "${params.outdir}/${sample_id}/nextclade/", pattern: "${sample_id}_${segment}_nextclade.json", mode: 'copy'
    publishDir "${params.outdir}/${sample_id}/nextclade/", pattern: "${sample_id}_${segment}_nextclade.tsv",  mode: 'copy'

    input:
    tuple val(sample_id), val(segment), path(consensus_fa), path(dataset_dir)

    output:
    tuple val(sample_id), val(segment), path("${sample_id}_${segment}_nextclade.json"), path(dataset_dir), emit: json_with_ref
    tuple val(sample_id), val(segment), path("${sample_id}_${segment}_nextclade.tsv"),                      emit: tsv

    script:
    """
    nextclade run \
        --input-dataset ${dataset_dir} \
        --output-json ${sample_id}_${segment}_nextclade.json \
        --output-tsv  ${sample_id}_${segment}_nextclade.tsv \
        ${consensus_fa}
    """
}


process nextclade_to_vcf {

    conda "${projectDir}/environments/nextclade.yml"

    tag { "${sample_id}_${segment}" }

    publishDir "${params.outdir}/${sample_id}/nextclade/", pattern: "${sample_id}_${segment}.vcf", mode: 'copy'

    input:
    tuple val(sample_id), val(segment), path(nextclade_json), path(dataset_dir)

    output:
    tuple val(sample_id), val(segment), path("${sample_id}_${segment}.vcf"), emit: vcf

    script:
    """
    nextclade_to_freebayes_vcf.py \
        -i ${nextclade_json} \
        -d ${dataset_dir} \
        -o ${sample_id}_${segment}.vcf \
        --sample ${sample_id}
    """
}
