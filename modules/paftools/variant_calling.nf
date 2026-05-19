// paftools/variant_calling.nf - paftools.js variant calling, sorting, annotation, and reformatting for the PAFTOOLS workflow

process paftools_call {

    tag { "${name}" }

    //publishDir "${params.outdir}/${sample_id}/vcf", pattern: "${name}_raw.vcf", mode: 'copy'

    input:
    tuple val(sample_id), val(segment), path(alignment_paf), path(reference_fasta)

    output:
    tuple val(sample_id), val(segment), path("${name}_raw.vcf"), env(SNPEFF_DB_NAME), emit: vcf
    tuple val(sample_id), path("${sample_id}_paftools_call_provenance.yml"), emit: provenance

    script:
    name = "${sample_id}_${segment}"
    """
    printf -- "- process_name: paftools_call\\n" > ${sample_id}_paftools_call_provenance.yml
    printf -- "  tools: \\n" >> ${sample_id}_paftools_call_provenance.yml
    printf -- "  - tool_name: paftools.js\\n    tool_version: \$(paftools.js version)\\n" >> ${sample_id}_paftools_call_provenance.yml

    paftools.js call -l 100 -L 100 -f ${reference_fasta} ${name}.paf > ${name}_raw.vcf

    SNPEFF_DB_NAME=\$(head -n1 ${reference_fasta} | cut -d'|' -f2)
    SNPEFF_DB_NAME="influenza_\${SNPEFF_DB_NAME,,}"

    echo "influenza_\${SNPEFF_DB_NAME}" > ${name}_snpeff_db_name.txt
    """
}


process paftools_sort_vcf {

    tag { name }

    //publishDir "${params.outdir}/${sample_id}/vcf", pattern: "${name}_sorted.vcf", mode: 'copy'

    input:
    tuple val(sample_id), val(segment), path(raw_vcf), val(snpeff_db_name)

    output:
    tuple val(sample_id), val(segment), path("${name}_sorted.vcf"), val(snpeff_db_name), emit: vcf
    tuple val(sample_id), path("${sample_id}_paftools_sort_vcf_provenance.yml"), emit: provenance

    script:
    name = "${sample_id}_${segment}"

    """
    printf -- "- process_name: paftools_sort_vcf\\n" > ${sample_id}_paftools_sort_vcf_provenance.yml
    printf -- "  tools: \\n  - tool_name: bcftools\\n    tool_version: \$(bcftools --version | head -n1 | cut -d' ' -f2)\\n" >> ${sample_id}_paftools_sort_vcf_provenance.yml

    bcftools sort ${raw_vcf} -o ${name}_sorted.vcf
    """
}


process paftools_snpeff_annotation {
    tag { name }

    conda "${projectDir}/environments/fluscan.yml"

    //publishDir "${params.outdir}/${sample_id}/vcf", pattern: "${sample_id}*__anno.vcf", mode: 'copy'

    input:
    tuple val(sample_id), val(segment), path(vcf), val(snpeff_db_name), path(snpeff_config_path)

    output:
    tuple val(sample_id), val(segment), path("${sample_id}*.vcf"), emit: vcf
    tuple val(sample_id), path("${sample_id}_paftools_snpeff_provenance.yml"), emit: provenance

    script:
    name = "${sample_id}_${segment}"

    """
    printf -- "- process_name: paftools_snpeff_annotation\\n" > ${sample_id}_paftools_snpeff_provenance.yml
    printf -- "  tools: \\n" >> ${sample_id}_paftools_snpeff_provenance.yml
    printf -- "  - tool_name: snpeff\\n    tool_version: \$(snpEff -version 2>&1 | cut -d\$'\\t' -f2)\\n" >> ${sample_id}_paftools_snpeff_provenance.yml
    printf -- "  - tool_name: bcftools\\n    tool_version: \$(bcftools --version | head -n1 | cut -d' ' -f2)\\n" >> ${sample_id}_paftools_snpeff_provenance.yml

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
