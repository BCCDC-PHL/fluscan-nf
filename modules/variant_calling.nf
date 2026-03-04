// variant_calling.nf - Variant calling and annotation processes

process freebayes {

    tag { sample_id }

    conda "${projectDir}/environments/fluscan.yml"

    publishDir "${params.outdir}/${sample_id}/vcf", pattern: "${sample_id}*{gz,csi}", mode: 'copy'

    input: 
	tuple val(sample_id), path(bamfile), path(bam_index), path(reference)
	
	output:
    tuple val(sample_id), path("${sample_id}*vcf.gz"), path("${sample_id}*.csi"), emit: vcf
    tuple val(sample_id), path("${sample_id}*provenance.yml"), emit: provenance

    script:

    """
	printf -- "- process_name: freebayes\\n"                                                                                           >> ${sample_id}_freebayes_provenance.yml
    printf -- "  tools: \\n  - tool_name: samtools\\n    tool_version: \$(samtools --version 2>&1 | head -n1 | cut -d' ' -f2)\\n"   >> ${sample_id}_freebayes_provenance.yml
    printf -- "  - tool_name: freebayes\\n    tool_version: \$(freebayes --version | cut -d' ' -f3)\\n"                             >> ${sample_id}_freebayes_provenance.yml
    printf -- "  - tool_name: bcftools\\n    tool_version: \$(bcftools --version | head -n1 | cut -d' ' -f2)\\n"                    >> ${sample_id}_freebayes_provenance.yml

    samtools faidx ${reference}
    freebayes -X -O \
    -b ${bamfile} \
    -v ${sample_id}.vcf \
    -f ${reference} \
    --limit-coverage 200 \
    -p 1 \
    -m 20 \
    -q 20 \
    --pooled-continuous \

    bcftools sort --output-type z ${sample_id}.vcf > ${sample_id}.vcf.gz
    bcftools index ${sample_id}.vcf.gz
    """

}


process snpeff_annotation {
    tag { sample_id }

    conda "${projectDir}/environments/fluscan.yml"

    publishDir "${params.outdir}/${sample_id}/vcf", pattern: "${sample_id}*{vcf,tsv}", mode: 'copy'

    input: 
	tuple val(sample_id), path(vcf), path(vcf_index), path(snpeff_config_path)
	
	output:
    tuple val(sample_id), path("${sample_id}*vcf"), emit: vcf
    tuple val(sample_id), path("${sample_id}*tsv"), emit: vcf_filter
    tuple val(sample_id), path("${sample_id}*provenance.yml"), emit: provenance

    script:

    """
	printf -- "- process_name: snpeff_annotation\\n"                                                      >> ${sample_id}_snpeff_provenance.yml
    printf -- "  tools: \\n  "                                                                            >> ${sample_id}_snpeff_provenance.yml
    printf -- "  - tool_name: snpeff\\n    tool_version: \$(snpEff -version 2>&1 | cut -d\$'\t' -f2)\\n"   >> ${sample_id}_snpeff_provenance.yml
    printf -- "  - tool_name: samtools\\n    tool_version: \$(samtools --version 2>&1 | head -n1 | cut -d' ' -f2)\\n"   >> ${sample_id}_snpeff_provenance.yml
    printf -- "  - tool_name: bcftools\\n    tool_version: \$(bcftools --version | head -n1 | cut -d' ' -f2)\\n"   >> ${sample_id}_snpeff_provenance.yml

    run_snpeff.py --input ${vcf} --config ${snpeff_config_path} --outname ${sample_id}

    reformat_vcf.py --input ${sample_id}_all_anno.vcf --output ${sample_id}_all_muts.tsv
    """

}


process mutation_watch {
    tag { sample_id }

    conda "${projectDir}/environments/fluscan.yml"

    publishDir "${params.outdir}/${sample_id}/vcf", pattern: "${sample_id}*{tsv}", mode: 'copy'

    input: 
	tuple val(sample_id), path(vcf_filter), path(watchlist)
	
	output:
    tuple val(sample_id), path("${sample_id}*tsv"), emit: vcf
    script:
    watchlist_name = watchlist.baseName
    """
    mutation_watcher.py --input ${vcf_filter} --watchlist ${watchlist} --output ${sample_id}_${watchlist_name}.tsv
    """

}
