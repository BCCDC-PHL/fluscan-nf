// variant_calling.nf - Variant calling and annotation processes

process lofreq {

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
	printf -- "- process_name: lofreq\\n"                                                                                              >> ${sample_id}_lofreq_provenance.yml
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


process consensus_paftools_call {

    tag { "${name}" }

    //publishDir "${params.outdir}/${sample_id}/vcf", pattern: "${name}_raw.vcf", mode: 'copy'

    input:
    tuple val(sample_id), val(segment), path(alignment_paf), path(reference_fasta) 

    output:
    tuple val(sample_id), val(segment), path("${name}_raw.vcf"), env(SNPEFF_DB_NAME), emit: vcf
    tuple val(sample_id), path("${sample_id}_consensus_paftools_provenance.yml"), emit: provenance

    script:
    name = "${sample_id}_${segment}"
    """
    printf -- "- process_name: consensus_paftools_call\\n" > ${sample_id}_consensus_paftools_provenance.yml
    printf -- "  tools: \\n" >> ${sample_id}_consensus_paftools_provenance.yml
    printf -- "  - tool_name: paftools.js\\n    tool_version: \$(paftools.js version)\\n" >> ${sample_id}_consensus_paftools_provenance.yml

    paftools.js call -l 100 -L 100 -f ${reference_fasta} ${name}.paf > ${name}_raw.vcf

    SNPEFF_DB_NAME=\$(head -n1 ${reference_fasta} | cut -d'|' -f2)
    SNPEFF_DB_NAME="influenza_\${SNPEFF_DB_NAME,,}"

    echo "influenza_\${SNPEFF_DB_NAME}" > ${name}_snpeff_db_name.txt
    """
}


process consensus_sort_vcf {

    tag { name }

    //publishDir "${params.outdir}/${sample_id}/vcf", pattern: "${name}_sorted.vcf", mode: 'copy'

    input:
    tuple val(sample_id), val(segment), path(raw_vcf), val(snpeff_db_name)

    output:
    tuple val(sample_id), val(segment), path("${name}_sorted.vcf"), val(snpeff_db_name), emit: vcf
    tuple val(sample_id), path("${sample_id}_consensus_sort_vcf_provenance.yml"), emit: provenance

    script:
    name = "${sample_id}_${segment}"

    """
    printf -- "- process_name: consensus_sort_vcf\\n" > ${sample_id}_consensus_sort_vcf_provenance.yml
    printf -- "  tools: \\n  - tool_name: bcftools\\n    tool_version: \$(bcftools --version | head -n1 | cut -d' ' -f2)\\n" >> ${sample_id}_consensus_sort_vcf_provenance.yml

    bcftools sort ${raw_vcf} -o ${name}_sorted.vcf
    """
}


process snpeff_annotation {
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
    printf -- "- process_name: snpeff_annotation\\n"                                                      >> ${sample_id}_snpeff_provenance.yml
    printf -- "  tools: \\n  "                                                                            >> ${sample_id}_snpeff_provenance.yml
    printf -- "  - tool_name: snpeff\\n    tool_version: \$(snpEff -version 2>&1 | cut -d\$'\t' -f2)\\n"   >> ${sample_id}_snpeff_provenance.yml
    printf -- "  - tool_name: samtools\\n    tool_version: \$(samtools --version 2>&1 | head -n1 | cut -d' ' -f2)\\n"   >> ${sample_id}_snpeff_provenance.yml
    printf -- "  - tool_name: bcftools\\n    tool_version: \$(bcftools --version | head -n1 | cut -d' ' -f2)\\n"   >> ${sample_id}_snpeff_provenance.yml

    run_snpeff.py --input ${vcf} --config ${snpeff_config_path} --outname ${sample_id}
    """

}


process consensus_snpeff_annotation {
    tag { name }

    conda "${projectDir}/environments/fluscan.yml"

    //publishDir "${params.outdir}/${sample_id}/vcf", pattern: "${sample_id}*__anno.vcf", mode: 'copy'

    input:
    tuple val(sample_id), val(segment), path(vcf), val(snpeff_db_name), path(snpeff_config_path)

    output:
    tuple val(sample_id), val(segment), path("${sample_id}*.vcf"), emit: vcf
    tuple val(sample_id), path("${sample_id}_consensus_snpeff_provenance.yml"), emit: provenance

    script:
    name = "${sample_id}_${segment}"
    
    """
    printf -- "- process_name: consensus_snpeff_annotation\\n" > ${sample_id}_consensus_snpeff_provenance.yml
    printf -- "  tools: \\n" >> ${sample_id}_consensus_snpeff_provenance.yml
    printf -- "  - tool_name: snpeff\\n    tool_version: \$(snpEff -version 2>&1 | cut -d\$'\\t' -f2)\\n" >> ${sample_id}_consensus_snpeff_provenance.yml
    printf -- "  - tool_name: bcftools\\n    tool_version: \$(bcftools --version | head -n1 | cut -d' ' -f2)\\n" >> ${sample_id}_consensus_snpeff_provenance.yml

    REF_NAME=\$(grep -oP "##contig=<ID=\\K([^,]+)" ${vcf})
    REF_NAME_FIELDS=(\${REF_NAME//|/ })

    awk 'BEGIN{OFS="\\t"} /^##contig=/{gsub(/\\|[^,>]*/, ""); print; next} /^#/{print; next} {split(\$1,a,"|"); \$1=a[1]; print}' ${vcf} \\
     > rename.vcf

    OUTPUT_FILENAME="${sample_id}__\${REF_NAME_FIELDS[0]}__\${REF_NAME_FIELDS[1]}__\${REF_NAME_FIELDS[2]}.vcf"

    snpEff -Xmx4g -noStats -hgvs1LetterAa -config ${snpeff_config_path} ${snpeff_db_name} rename.vcf > \${OUTPUT_FILENAME}
    """
}

process paftools_to_freebayes_vcf {
    tag { "${sample_id}_${segment}" }

    //publishDir "${params.outdir}/${sample_id}/vcf", pattern: "${sample_id}*_reformat.vcf", mode: 'copy'

    input:
    tuple val(sample_id), val(segment), path(paftools_vcf), path(coverage_tsv)

    output:
    tuple val(sample_id), val(segment), path("${sample_id}*.vcf"), emit: vcf

    script:
    """
    name=\$(basename ${paftools_vcf} .vcf)

    paftools_to_freebayes_vcf.py --input ${paftools_vcf} --output tmp.vcf

    # Inject alignment-coverage status into the VCF header
    COV_STATUS=\$(awk 'NR==2 {print \$11}' ${coverage_tsv})
    COV_FRAC=\$(awk 'NR==2 {print \$6}' ${coverage_tsv})
    COV_FIRST=\$(awk 'NR==2 {print \$7}' ${coverage_tsv})
    COV_LAST=\$(awk 'NR==2 {print \$8}' ${coverage_tsv})

    awk -v hdr="##fluscan_alignment_coverage=STATUS=\${COV_STATUS};fraction=\${COV_FRAC};range=\${COV_FIRST}-\${COV_LAST}" \\
      'FNR==1 && /^##fileformat/{print; print hdr; next} {print}' \\
      tmp.vcf > \${name}_paf.vcf
    """
}


process post_process_vcf {
    tag { "${sample_id}_${segment}" }

    conda "${projectDir}/environments/fluscan.yml"

    publishDir "${params.outdir}/${sample_id}/${workflow}", pattern: "${sample_id}*_final.vcf", mode: 'copy'

    input:
    tuple val(sample_id), val(segment), path(vcf), path(reference)

    output:
    tuple val(sample_id), val(segment), path("${sample_id}*_final.vcf"), emit: vcf
    tuple val(sample_id), val(segment), path("${sample_id}_post_process_vcf_provenance.yml"), emit: provenance

    script:
    def min_af = params.mnv_min_af ?: 0.5
    workflow = task.ext.workflow ?: "unknown_workflow"
    """
    printf -- "- process_name: post_process_vcf\\n" > ${sample_id}_post_process_vcf_provenance.yml
    printf -- "  tools: \\n" >> ${sample_id}_post_process_vcf_provenance.yml
    printf -- "  - tool_name: python\\n    tool_version: \$(python --version 2>&1 | cut -d' ' -f2)\\n" >> ${sample_id}_post_process_vcf_provenance.yml
    printf -- "  - tool_name: fix_mnv_codons.py\\n    tool_version: 1.0.0\\n" >> ${sample_id}_post_process_vcf_provenance.yml

    name=\$(basename ${vcf} .vcf)

    echo "Processing ${vcf} with min_af=${min_af}"

    echo "Renaming reference..."

    sed -i -E 's/^>([^\\|]+)\\|.+/>\\1/g' ${reference}

    bcftools norm -f ${reference} -o ${name}_norm.vcf ${vcf}

    fix_mnv_codons.py \\
    --input ${name}_norm.vcf \\
    --reference ${reference} \\
    --output-vcf \${name}_final.vcf \\
    --min-af ${min_af}
    """
}



process summarize_vcf {
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
    printf -- "- process_name: summarize_vcf\\n" > ${sample_id}_summarize_vcf_provenance.yml
    printf -- "  tools: \\n" >> ${sample_id}_summarize_vcf_provenance.yml
    printf -- "  - tool_name: python\\n    tool_version: \$(python --version 2>&1 | cut -d' ' -f2)\\n" >> ${sample_id}_summarize_vcf_provenance.yml
    printf -- "  - tool_name: summarize_vcf.py\\n    tool_version: 1.0.0\\n" >> ${sample_id}_summarize_vcf_provenance.yml

    summarize_vcf.py --input ${vcf} --output ${sample_id}_mutations.tsv
    """
}


process mutation_watch {
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
