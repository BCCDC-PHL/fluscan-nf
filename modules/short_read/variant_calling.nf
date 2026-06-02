// short_read/variant_calling.nf - Variant calling and annotation processes for the SHORT_READ workflow

process short_read_lofreq {

    tag { sample_id }

    conda "${projectDir}/environments/fluscan.yml"

    //publishDir "${params.outdir}/${sample_id}/vcf", pattern: "${sample_id}*{gz,csi}", mode: 'copy'

    input:
	tuple val(sample_id), path(bamfile), path(bam_index), path(reference)

	output:
    tuple val(sample_id), path("${sample_id}*vcf.gz"), path("${sample_id}*.csi"), emit: vcf
    tuple val(sample_id), path("${sample_id}*provenance.yml"), emit: provenance

    script:

    """
	printf -- "- process_name: short_read_lofreq\\n"                                                                                           >> ${sample_id}_lofreq_provenance.yml
    printf -- "  tools: \\n  - tool_name: samtools\\n    tool_version: \$(samtools --version 2>&1 | head -n1 | cut -d' ' -f2)\\n"   >> ${sample_id}_lofreq_provenance.yml
    printf -- "  - tool_name: lofreq\\n    tool_version: \$(lofreq version 2>&1 | head -n1)\\n"                                     >> ${sample_id}_lofreq_provenance.yml
    printf -- "  - tool_name: bcftools\\n    tool_version: \$(bcftools --version | head -n1 | cut -d' ' -f2)\\n"                    >> ${sample_id}_lofreq_provenance.yml

    samtools faidx ${reference}

    lofreq viterbi -f ${reference} ${bamfile} | samtools sort -o ${sample_id}_viterbi.bam
    samtools index ${sample_id}_viterbi.bam

    lofreq indelqual --dindel -f ${reference} ${sample_id}_viterbi.bam -o ${sample_id}_indelqual.bam
    samtools index ${sample_id}_indelqual.bam

    TARGET_COV=${params.lofreq_target_cov}

    MEAN_DEPTH=\$(samtools coverage ${sample_id}_indelqual.bam | awk 'NR>1 {sum+=\$7; count++} END {if(count>0) print sum/count; else print 0}')

    FRACTION=\$(echo "scale=4; \$TARGET_COV / \$MEAN_DEPTH" | bc)

    echo "Target coverage: \$TARGET_COV"
    echo "Mean depth: \$MEAN_DEPTH"
    echo "Subsampling fraction: \$FRACTION"

    if (( \$(echo "\$FRACTION < 1" | bc -l) )); then
        echo "Subsampling to target coverage..."
        samtools view -s "0\${FRACTION}" -b ${sample_id}_indelqual.bam > ${sample_id}_final.bam
    else
        echo "Mean depth is below target coverage; skipping subsampling."
        cp ${sample_id}_indelqual.bam ${sample_id}_final.bam
    fi

    samtools index ${sample_id}_final.bam

    lofreq call-parallel --pp-threads 8 -f ${reference} -o ${sample_id}.vcf ${sample_id}_final.bam

    bcftools reheader --fai ${reference}.fai ${sample_id}.vcf | bcftools sort --output-type z > ${sample_id}.vcf.gz
    bcftools index ${sample_id}.vcf.gz
    """

}


process short_read_snpeff_annotation {
    tag { sample_id }

    conda "${projectDir}/environments/fluscan.yml"

    //publishDir "${params.outdir}/${sample_id}/vcf", pattern: "${sample_id}*__anno.vcf", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf), path(vcf_index), path(snpeff_config_path)

    output:
    tuple val(sample_id), path("${sample_id}_all_anno.vcf"), emit: all_vcf
    tuple val(sample_id), path("${sample_id}*__anno.vcf"), emit: vcfs
    tuple val(sample_id), path("${sample_id}_snpeff_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: short_read_snpeff_annotation\\n"                                                 >> ${sample_id}_snpeff_provenance.yml
    printf -- "  tools: \\n  "                                                                            >> ${sample_id}_snpeff_provenance.yml
    printf -- "  - tool_name: snpeff\\n    tool_version: \$(snpEff -version 2>&1 | cut -d\$'\t' -f2)\\n"   >> ${sample_id}_snpeff_provenance.yml
    printf -- "  - tool_name: samtools\\n    tool_version: \$(samtools --version 2>&1 | head -n1 | cut -d' ' -f2)\\n"   >> ${sample_id}_snpeff_provenance.yml
    printf -- "  - tool_name: bcftools\\n    tool_version: \$(bcftools --version | head -n1 | cut -d' ' -f2)\\n"   >> ${sample_id}_snpeff_provenance.yml

    run_snpeff.py --input ${vcf} --config ${snpeff_config_path} --outname ${sample_id}
    """

}


process short_read_post_process_vcf {
    tag { "${sample_id}:${segment}" }

    conda "${projectDir}/environments/fluscan.yml"

    publishDir "${params.outdir}/${sample_id}/short_read", pattern: "${sample_id}*_final.vcf", mode: 'copy'

    input:
    tuple val(sample_id), val(segment), path(vcf), path(reference)

    output:
    tuple val(sample_id), path("${sample_id}*_final.vcf"), emit: vcf
    tuple val(sample_id), path("${sample_id}_${segment}_post_process_vcf_provenance.yml"), emit: provenance

    script:
    def min_af = params.mnv_min_af ?: 0.5
    """
    printf -- "- process_name: short_read_post_process_vcf\\n" > ${sample_id}_${segment}_post_process_vcf_provenance.yml
    printf -- "  tools: \\n" >> ${sample_id}_${segment}_post_process_vcf_provenance.yml
    printf -- "  - tool_name: python\\n    tool_version: \$(python --version 2>&1 | cut -d' ' -f2)\\n" >> ${sample_id}_${segment}_post_process_vcf_provenance.yml
    printf -- "  - tool_name: fix_mnv_codons.py\\n    tool_version: 1.0.0\\n" >> ${sample_id}_${segment}_post_process_vcf_provenance.yml

    sed -E 's/^>([^|]+)\\|.+/>\\1/g' ${reference} > reference_cleaned.fasta

    vcf_basename=\$(basename ${vcf} .vcf)

    fix_mnv_codons.py \
        --input ${vcf} \
        --reference reference_cleaned.fasta \
        --output-vcf \${vcf_basename}_final.vcf \
        --min-af ${min_af}
    """
}


process short_read_summarize_vcf {
    tag { sample_id }

    errorStrategy 'ignore'

    conda "${projectDir}/environments/fluscan.yml"

    // publishDir "${params.outdir}/${sample_id}/vcf", pattern: "${sample_id}*tsv", mode: 'copy'

    input:
    tuple val(sample_id), path(vcf)

    output:
    tuple val(sample_id), path("${sample_id}_mutations.tsv"), emit: tsv
    tuple val(sample_id), path("${sample_id}_summarize_vcf_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: short_read_summarize_vcf\\n" > ${sample_id}_summarize_vcf_provenance.yml
    printf -- "  tools: \\n" >> ${sample_id}_summarize_vcf_provenance.yml
    printf -- "  - tool_name: python\\n    tool_version: \$(python --version 2>&1 | cut -d' ' -f2)\\n" >> ${sample_id}_summarize_vcf_provenance.yml
    printf -- "  - tool_name: summarize_vcf.py\\n    tool_version: 1.0.0\\n" >> ${sample_id}_summarize_vcf_provenance.yml

    summarize_vcf.py --input ${vcf} --output ${sample_id}_mutations.tsv
    """
}


process short_read_mutation_watch {
    tag { sample_id }

    conda "${projectDir}/environments/fluscan.yml"

    //publishDir "${params.outdir}/${sample_id}/vcf", pattern: "${sample_id}*{tsv}", mode: 'copy'

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
