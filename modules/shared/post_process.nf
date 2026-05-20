// shared/post_process.nf - MNV codon correction applied to annotated VCFs (shared by all workflows)

process post_process_vcf {
    tag { "${sample_id}_${segment}" }

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

    fix_mnv_codons.py \\
    --input ${vcf} \\
    --reference ${reference} \\
    --output-vcf \${name}_final.vcf \\
    --min-af ${min_af}
    """
}
