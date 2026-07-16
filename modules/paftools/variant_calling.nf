// paftools/variant_calling.nf - paftools.js variant calling, sorting, annotation, and reformatting for the PAFTOOLS workflow

process paftools_call {

    tag { sample_id }

    input:
    tuple val(sample_id),
          path(pafs, stageAs: 'pafs/*'),
          path(references, stageAs: 'ref_inputs/*')

    output:
    tuple val(sample_id), path('*_raw.vcf'), path('reference/*'), emit: vcf
    tuple val(sample_id), path("${sample_id}_paftools_call_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: paftools_call\\n" > ${sample_id}_paftools_call_provenance.yml
    printf -- "  tools: \\n" >> ${sample_id}_paftools_call_provenance.yml
    printf -- "  - tool_name: paftools.js\\n    tool_version: \$(paftools.js version)\\n" >> ${sample_id}_paftools_call_provenance.yml

    for paf_file in pafs/*.paf; do
        paf_basename=\$(basename "\$paf_file" .paf)
        segment="\${paf_basename#${sample_id}_}"
        reference_fa="ref_inputs/${sample_id}_ref_\${segment}.s.fasta"
        output_prefix="${sample_id}_\${segment}"
        paftools.js call -l 100 -L 100 -f "\$reference_fa" "\$paf_file" > "\${output_prefix}_raw.vcf"
    done
    cp -r ref_inputs reference
    """
}


process paftools_sort_vcf {

    tag { sample_id }

    input:
    tuple val(sample_id),
          path(raw_vcfs, stageAs: 'vcf/*'),
          path(references, stageAs: 'ref_inputs/*')

    output:
    tuple val(sample_id), path('*_sorted.vcf'), path('reference/*'), emit: vcf
    tuple val(sample_id), path("${sample_id}_paftools_sort_vcf_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: paftools_sort_vcf\\n" > ${sample_id}_paftools_sort_vcf_provenance.yml
    printf -- "  tools: \\n  - tool_name: bcftools\\n    tool_version: \$(bcftools --version | head -n1 | cut -d' ' -f2)\\n" >> ${sample_id}_paftools_sort_vcf_provenance.yml

    for vcf_file in vcf/*_raw.vcf; do
        output_prefix=\$(basename "\$vcf_file" _raw.vcf)
        bcftools sort "\$vcf_file" -o "\${output_prefix}_sorted.vcf"
    done
    cp -r ref_inputs reference
    """
}


process paftools_snpeff_annotation {
    tag { sample_id }

    conda "${projectDir}/environments/fluscan.yml"

    input:
    tuple val(sample_id),
          path(sorted_vcfs, stageAs: 'vcf/*'),
          path(references, stageAs: 'reference/*'),
          path(snpeff_config_path)

    output:
    tuple val(sample_id), path("${sample_id}*.vcf"), emit: vcf
    tuple val(sample_id), path("${sample_id}_paftools_snpeff_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: paftools_snpeff_annotation\\n" > ${sample_id}_paftools_snpeff_provenance.yml
    printf -- "  tools: \\n" >> ${sample_id}_paftools_snpeff_provenance.yml
    printf -- "  - tool_name: snpeff\\n    tool_version: \$(snpEff -version 2>&1 | cut -d\$'\\t' -f2)\\n" >> ${sample_id}_paftools_snpeff_provenance.yml
    printf -- "  - tool_name: bcftools\\n    tool_version: \$(bcftools --version | head -n1 | cut -d' ' -f2)\\n" >> ${sample_id}_paftools_snpeff_provenance.yml

    for vcf_file in vcf/*_sorted.vcf; do
        vcf_basename=\$(basename "\$vcf_file" _sorted.vcf)
        segment="\${vcf_basename#${sample_id}_}"
        reference_fa="reference/${sample_id}_ref_\${segment}.s.fasta"

        REF_NAME=\$(grep -oP "##contig=<ID=\\K([^,]+)" "\$vcf_file")
        REF_NAME_FIELDS=(\${REF_NAME//|/ })

        SNPEFF_DB_NAME=\$(head -n1 "\$reference_fa" | cut -d'|' -f2)
        SNPEFF_DB_NAME="influenza_\${SNPEFF_DB_NAME,,}"

        awk 'BEGIN{OFS="\\t"} /^##contig=/{gsub(/\\|[^,>]*/, ""); print; next} /^#/{print; next} {split(\$1,a,"|"); \$1=a[1]; print}' "\$vcf_file" \\
         > rename_\${segment}.vcf

        OUTPUT_FILENAME="${sample_id}__\${REF_NAME_FIELDS[0]}__\${REF_NAME_FIELDS[1]}__\${REF_NAME_FIELDS[2]}.vcf"

        snpEff -Xmx4g -noStats -hgvs1LetterAa -config ${snpeff_config_path} "\$SNPEFF_DB_NAME" rename_\${segment}.vcf > "\${OUTPUT_FILENAME}"
    done
    """
}


process paftools_to_freebayes_vcf {
    tag { sample_id }

    input:
    tuple val(sample_id),
          path(anno_vcfs, stageAs: 'vcf/*'),
          path(coverage_tsvs, stageAs: 'cov/*')

    output:
    tuple val(sample_id), path("${sample_id}*_paf.vcf"), emit: vcf

    script:
    """
    for anno_vcf in vcf/*.vcf; do
        vcf_basename=\$(basename "\$anno_vcf" .vcf)
        segment=\$(echo "\$vcf_basename" | awk -F'__' '{print \$NF}')
        cov_tsv="cov/${sample_id}_\${segment}_coverage_qc.tsv"

        paftools_to_freebayes_vcf.py --input "\$anno_vcf" --output tmp_\${segment}.vcf

        COV_STATUS=\$(awk 'NR==2 {print \$11}' "\$cov_tsv")
        COV_FRAC=\$(awk 'NR==2 {print \$6}' "\$cov_tsv")
        COV_FIRST=\$(awk 'NR==2 {print \$7}' "\$cov_tsv")
        COV_LAST=\$(awk 'NR==2 {print \$8}' "\$cov_tsv")

        awk -v coverage_header="##fluscan_alignment_coverage=STATUS=\${COV_STATUS};fraction=\${COV_FRAC};range=\${COV_FIRST}-\${COV_LAST}" \\
          'FNR==1 && /^##fileformat/{print; print coverage_header; next} {print}' \\
          tmp_\${segment}.vcf > \${vcf_basename}_paf.vcf
    done
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
