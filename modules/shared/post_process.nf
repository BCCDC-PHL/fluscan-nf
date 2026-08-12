// shared/post_process.nf - MNV codon correction applied to annotated VCFs (shared by all workflows)

// Per-sample batched form: processes all 8 segments in one task (paftools workflow).
// VCF files named {sample_id}__{ref}__{subtype}__{segment}_paf.vcf;
// references named {sample_id}_ref_{segment}.s.fasta — matched by segment in the loop.
process post_process_vcf {
    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}/${workflow}", pattern: "${sample_id}*_final.vcf", mode: 'copy'

    input:
    tuple val(sample_id),
          path(vcfs, stageAs: 'vcf/*'),
          path(references, stageAs: 'reference/*')

    output:
    tuple val(sample_id), path("${sample_id}*_final.vcf"), emit: vcf
    tuple val(sample_id), path("${sample_id}_post_process_vcf_provenance.yml"), emit: provenance

    script:
    def min_af = params.mnv_min_af ?: 0.5
    def parallel_jobs = task.cpus ?: 1
    workflow = task.ext.workflow ?: "unknown_workflow"
    """
    printf -- "- process_name: post_process_vcf\\n" > ${sample_id}_post_process_vcf_provenance.yml
    printf -- "  tools: \\n" >> ${sample_id}_post_process_vcf_provenance.yml
    printf -- "  - tool_name: python\\n    tool_version: \$(python --version 2>&1 | cut -d' ' -f2)\\n" >> ${sample_id}_post_process_vcf_provenance.yml
    printf -- "  - tool_name: parallel\\n    tool_version: \$(parallel --version | head -n1 | cut -d' ' -f3)\\n" >> ${sample_id}_post_process_vcf_provenance.yml
    printf -- "  - tool_name: fix_mnv_codons.py\\n    tool_version: 1.0.0\\n" >> ${sample_id}_post_process_vcf_provenance.yml

    for reference_fa in reference/*.fasta; do
        sed -i -E 's/^>([^|]+)\\|.+/>\\1/g' "\$reference_fa"
    done

    cat > run_fix_mnv_sample.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

vcf_file="\$1"
sample_id="\$2"
min_af="\$3"

vcf_basename=\$(basename "\$vcf_file" .vcf)
segment="\${vcf_basename##*__}"
segment="\${segment%_paf}"
reference_fa="reference/\${sample_id}_ref_\${segment}.s.fasta"

echo "Processing \$vcf_file with min_af=\${min_af}"

fix_mnv_codons.py \
    --input "\$vcf_file" \
    --reference "\$reference_fa" \
    --output-vcf "\${vcf_basename}_final.vcf" \
    --min-af "\${min_af}"
EOF
    chmod +x run_fix_mnv_sample.sh

    find vcf -maxdepth 1 -name '*.vcf' -print0 | \
        parallel -0 --will-cite --jobs ${parallel_jobs} --halt soon,fail=1 \
        ./run_fix_mnv_sample.sh {} '${sample_id}' '${min_af}'
    """
}


// Per-segment cross-sample batch form: processes all samples for one segment in one task (nextclade workflow).
// VCF files named {sample_id}__{ref}__{subtype}__{segment}_nc.vcf; all share the segment reference.
// publishDir uses saveAs to route each sample's VCF to its own output directory.
process post_process_vcf_batch {
    tag { segment }

    publishDir "${params.outdir}", mode: 'copy',
        saveAs: { filename ->
            def workflow  = task.ext.workflow ?: 'unknown_workflow'
            def sample_id = filename.split('__')[0]
            // Checked first: a no-mature VCF also ends in _final.vcf, and belongs in its
            // own subfolder rather than beside the VCF the rest of the workflow consumes.
            if (filename.endsWith('_nomat_final.vcf')) {
                return "${sample_id}/${workflow}/no_mature_peptide/${filename}"
            }
            if (filename.endsWith('_final.vcf')) {
                return "${sample_id}/${workflow}/${filename}"
            }
            return null
        }

    input:
    tuple val(segment),
          path(vcfs, stageAs: 'vcf/*'),
          path(dataset_dir)

    output:
    tuple val(segment), path("*_final.vcf"), emit: vcf
    tuple val(segment), path("post_process_${segment}_provenance.yml"), emit: provenance

    script:
    def min_af = params.mnv_min_af ?: 0.5
    def parallel_jobs = task.cpus ?: 1
    """
    printf -- "- process_name: post_process_vcf\\n" > post_process_${segment}_provenance.yml
    printf -- "  tools: \\n" >> post_process_${segment}_provenance.yml
    printf -- "  - tool_name: python\\n    tool_version: \$(python --version 2>&1 | cut -d' ' -f2)\\n" >> post_process_${segment}_provenance.yml
    printf -- "  - tool_name: parallel\\n    tool_version: \$(parallel --version | head -n1 | cut -d' ' -f3)\\n" >> post_process_${segment}_provenance.yml
    printf -- "  - tool_name: fix_mnv_codons.py\\n    tool_version: 1.0.0\\n" >> post_process_${segment}_provenance.yml

    cp ${dataset_dir}/reference.fasta reference_renamed.fasta
    sed -i -E 's/^>([^|]+)\\|.+/>\\1/g' reference_renamed.fasta

    cat > run_fix_mnv_batch.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

vcf_file="\$1"
min_af="\$2"

vcf_basename=\$(basename "\$vcf_file" .vcf)

echo "Processing \$vcf_file with min_af=\${min_af}"

fix_mnv_codons.py \
    --input "\$vcf_file" \
    --reference reference_renamed.fasta \
    --output-vcf "\${vcf_basename}_final.vcf" \
    --min-af "\${min_af}"
EOF
    chmod +x run_fix_mnv_batch.sh

    find vcf -maxdepth 1 -name '*.vcf' -print0 | \
        parallel -0 --will-cite --jobs ${parallel_jobs} --halt soon,fail=1 \
        ./run_fix_mnv_batch.sh {} '${min_af}'
    """
}
