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


process consensus_paftools_call {

    tag { sample_id }

    conda "${projectDir}/environments/fluscan.yml"

    publishDir "${params.outdir}/${sample_id}/vcf", pattern: "${sample_id}_raw.vcf", mode: 'copy'

    input:
    tuple val(sample_id), path(consensus_fasta), path(reference_fasta), val(subtype), val(snpeff_db_name)

    output:
    tuple val(sample_id), path("${sample_id}_raw.vcf"), val(snpeff_db_name), emit: vcf
    tuple val(sample_id), path("${sample_id}_consensus_paftools_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: consensus_paftools_call\\n" > ${sample_id}_consensus_paftools_provenance.yml
    printf -- "  tools: \\n" >> ${sample_id}_consensus_paftools_provenance.yml
    printf -- "  - tool_name: minimap2\\n    tool_version: \$(minimap2 --version 2>&1 | head -n1)\\n" >> ${sample_id}_consensus_paftools_provenance.yml
    printf -- "  - tool_name: k8\\n    tool_version: \$(k8 -v 2>&1 | head -n1)\\n" >> ${sample_id}_consensus_paftools_provenance.yml

    paftools_js=\$(command -v paftools.js || true)
    if [ -z "\$paftools_js" ]; then
        paftools_js=\$(find "\$(dirname "\$(command -v minimap2)")" -maxdepth 3 -name paftools.js | head -n1)
    fi
    if [ -z "\$paftools_js" ]; then
        echo "ERROR: Unable to locate paftools.js in the active environment." >&2
        exit 1
    fi

    printf -- "  - tool_name: paftools.js\\n    tool_version: \$(basename "\$paftools_js")\\n" >> ${sample_id}_consensus_paftools_provenance.yml

    minimap2 -cx asm5 --cs ${reference_fasta} ${consensus_fasta} | \
        k8 "\$paftools_js" call -l 1 -L 1 -f ${reference_fasta} - \
        > ${sample_id}_raw.vcf
    """
}


process consensus_sort_vcf {

    tag { sample_id }

    conda "${projectDir}/environments/fluscan.yml"

    publishDir "${params.outdir}/${sample_id}/vcf", pattern: "${sample_id}_sorted.vcf", mode: 'copy'

    input:
    tuple val(sample_id), path(raw_vcf), val(snpeff_db_name)

    output:
    tuple val(sample_id), path("${sample_id}_sorted.vcf"), val(snpeff_db_name), emit: vcf
    tuple val(sample_id), path("${sample_id}_consensus_sort_vcf_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: consensus_sort_vcf\\n" > ${sample_id}_consensus_sort_vcf_provenance.yml
    printf -- "  tools: \\n  - tool_name: bcftools\\n    tool_version: \$(bcftools --version | head -n1 | cut -d' ' -f2)\\n" >> ${sample_id}_consensus_sort_vcf_provenance.yml

    bcftools sort ${raw_vcf} -o ${sample_id}_sorted.vcf
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


process consensus_snpeff_annotation {
    tag { sample_id }

    conda "${projectDir}/environments/fluscan.yml"

    publishDir "${params.outdir}/${sample_id}/vcf", pattern: "${sample_id}_annotated.vcf", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf), val(snpeff_db_name), path(snpeff_config_path)

    output:
    tuple val(sample_id), path("${sample_id}_annotated.vcf"), emit: vcf
    tuple val(sample_id), path("${sample_id}_consensus_snpeff_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: consensus_snpeff_annotation\\n" > ${sample_id}_consensus_snpeff_provenance.yml
    printf -- "  tools: \\n" >> ${sample_id}_consensus_snpeff_provenance.yml
    printf -- "  - tool_name: snpeff\\n    tool_version: \$(snpEff -version 2>&1 | cut -d\$'\\t' -f2)\\n" >> ${sample_id}_consensus_snpeff_provenance.yml
    printf -- "  - tool_name: bcftools\\n    tool_version: \$(bcftools --version | head -n1 | cut -d' ' -f2)\\n" >> ${sample_id}_consensus_snpeff_provenance.yml

    snpEff -Xmx4g -noStats -config ${snpeff_config_path} ${snpeff_db_name} ${vcf} > ${sample_id}_annotated.vcf
    """
}


process consensus_reformat_vcf {
    tag { sample_id }

    conda "${projectDir}/environments/fluscan.yml"

    publishDir "${params.outdir}/${sample_id}/vcf", pattern: "${sample_id}_final_summary.tsv", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf)

    output:
    tuple val(sample_id), path("${sample_id}_final_summary.tsv"), emit: summary
    tuple val(sample_id), path("${sample_id}_consensus_reformat_vcf_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: consensus_reformat_vcf\\n" > ${sample_id}_consensus_reformat_vcf_provenance.yml
    printf -- "  tools: \\n" >> ${sample_id}_consensus_reformat_vcf_provenance.yml
    printf -- "  - tool_name: python\\n    tool_version: \$(python --version 2>&1 | cut -d' ' -f2)\\n" >> ${sample_id}_consensus_reformat_vcf_provenance.yml
    printf -- "  - tool_name: reformat_vcf.py\\n    tool_version: 1.0.0\\n" >> ${sample_id}_consensus_reformat_vcf_provenance.yml

    reformat_vcf.py --input ${vcf} --output ${sample_id}_final_summary.tsv
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
