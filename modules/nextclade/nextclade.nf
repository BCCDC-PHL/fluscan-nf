// nextclade/nextclade.nf - Per-segment Nextclade alignment and VCF conversion processes

// Cross-sample per-segment batch: runs nextclade once for one bundle's samples' sequences
// for a given segment and converts it to one final VCF per sample, with one QC row per sample.
// A sample that fails QC, including one Nextclade cannot align, gets no VCF.
process nextclade_run_batch {

    conda "${projectDir}/environments/nextclade.yml"

    tag { "${db_name}/${segment}" }

    publishDir "${params.outdir}", mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith('_final.vcf')) {
                def sample_id = filename.split('__')[0]
                def subdir    = filename.endsWith('_nc_nomat_final.vcf') ? 'no_mature_peptide/' : ''
                return "${sample_id}/nextclade/${subdir}${filename}"
            }
            // Not named ${segment} in the work dir: the staged dataset already holds that name.
            return filename == 'nc_output' ? "aggregate/nextclade/${db_name}/${segment}" : null
        }

    input:
    tuple val(db_name), val(segment), path(batch_fasta, stageAs: 'batch.fasta'), path(dataset_dir)

    output:
    tuple val(db_name), val(segment), path("*_nc_final.vcf"),              optional: true, emit: vcf
    tuple val(db_name), val(segment), path("*_nc_nomat_final.vcf"),        optional: true, emit: vcf_no_mature
    tuple val(db_name), val(segment), path("${segment}_nextclade_qc.tsv"),                 emit: qc
    tuple val(db_name), val(segment), path("nc_output")
    tuple val(db_name), val(segment), path("${segment}_nextclade_provenance.yml"),         emit: provenance

    script:
    def threads = task.cpus ?: 1
    def stem = "__\${REFERENCE_NAME}__\${SUBTYPE}__${segment}"
    def no_mature_arg = params.emit_no_mature_vcf ? "--no-mature-suffix \"${stem}_nc_nomat_final.vcf\"" : ""
    """
    nextclade run \\
        --input-dataset ${dataset_dir} \\
        --output-all nc_output \\
        --output-basename batch \\
        --include-reference \\
        --alignment-preset ${params.nextclade_alignment_preset} \\
        --jobs ${threads} \\
        batch.fasta

    REFERENCE_NAME=\$(head -n1 ${dataset_dir}/reference.fasta | sed -E 's/>(\\S+).*/\\1/g')
    SUBTYPE=\$(head -n1 ${dataset_dir}/reference.fasta | sed -E 's/.+\\(\\S+\\((\\S+)\\)\\S*\\).+/\\1/g')
    if echo "\$SUBTYPE" | grep -q '^>'; then SUBTYPE="unknown"; fi

    printf -- "- process_name: nextclade_run_batch\\n"                                           > ${segment}_nextclade_provenance.yml
    printf -- "  tools:\\n"                                                                       >> ${segment}_nextclade_provenance.yml
    printf -- "    - tool_name: python\\n      tool_version: \$(python3 --version 2>&1 | cut -d' ' -f2)\\n" >> ${segment}_nextclade_provenance.yml
    printf -- "    - tool_name: nextclade\\n      tool_version: \$(nextclade --version | awk '{print \$NF}')\\n" >> ${segment}_nextclade_provenance.yml
    printf -- "  parameters:\\n"                                                                  >> ${segment}_nextclade_provenance.yml
    printf -- "    - parameter: segment\\n      value: ${segment}\\n"                             >> ${segment}_nextclade_provenance.yml
    printf -- "    - parameter: reference\\n      value: \$REFERENCE_NAME\\n"                     >> ${segment}_nextclade_provenance.yml
    printf -- "    - parameter: alignment_preset\\n      value: ${params.nextclade_alignment_preset}\\n" >> ${segment}_nextclade_provenance.yml
    printf -- "    - parameter: consensus_coverage_pass_pct\\n      value: ${params.consensus_coverage_pass_pct}\\n" >> ${segment}_nextclade_provenance.yml

    nextclade_to_freebayes_vcf.py \\
        -i nc_output/batch.json \\
        -d ${dataset_dir} \\
        --suffix "${stem}_nc_final.vcf" \\
        ${no_mature_arg} \\
        --segment ${segment} \\
        --qc ${segment}_nextclade_qc.tsv \\
        --pass-pct ${params.consensus_coverage_pass_pct}
    """
}
