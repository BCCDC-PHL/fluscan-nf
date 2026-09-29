// nextclade/nextclade.nf - Per-segment Nextclade alignment and VCF conversion processes

// Cross-sample per-segment batch: runs nextclade once for all samples' sequences for a
// given segment, splits the output into per-sample JSON and TSV files, and converts it to
// one final VCF per sample. A sample Nextclade cannot align gets no JSON and no VCF.
process nextclade_run_batch {

    conda "${projectDir}/environments/nextclade.yml"

    tag { segment }

    publishDir "${params.outdir}", mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith('_final.vcf')) {
                def workflow  = task.ext.workflow ?: 'unknown_workflow'
                def sample_id = filename.split('__')[0]
                def subdir    = filename.endsWith('_nc_nomat_final.vcf') ? 'no_mature_peptide/' : ''
                return "${sample_id}/${workflow}/${subdir}${filename}"
            }
            if (filename =~ /.*_nextclade_status\.tsv$/) {
                def workflow = task.ext.workflow ?: 'unknown_workflow'
                return "${workflow}/status/${filename}"
            }
            if (filename =~ /.*_nextclade\.(json|tsv)$/) {
                def file_extension = filename.endsWith('.json') ? 'json' : 'tsv'
                def filename_base  = filename.replaceAll(/_nextclade\.(json|tsv)$/, '')
                def sample_id      = (filename =~ /^[^_]+/)[0]
                def workflow       = task.ext.workflow ?: 'unknown_workflow'
                return "${sample_id}/${workflow}/${file_extension}/${filename}"
            }
            if (filename =~ /nc_output.*/) {
                def newname = filename.replaceAll(/nc_output/, "${segment}_output")
                return "nextclade/${newname}"
            }
            return null
        }

    input:
    tuple val(segment), val(sample_ids),
          path(consensus_fastas, stageAs: 'inputs/*'),
          path(dataset_dir)

    output:
    tuple val(segment), path("*_nc_final.vcf"),              optional: true, emit: vcf
    tuple val(segment), path("*_nc_nomat_final.vcf"),        optional: true, emit: vcf_no_mature
    tuple val(segment), path("*_nextclade.json"),            optional: true, emit: json
    tuple val(segment), path("*_nextclade.tsv"),                             emit: tsv
    tuple val(segment), path("nc_output"),                                   emit: all
    tuple val(segment), path("${segment}_nextclade_status.tsv"),             emit: status

    script:
    def threads = task.cpus ?: 1
    def stem = "__\${REFERENCE_NAME}__\${SUBTYPE}__${segment}"
    def no_mature_arg = params.emit_no_mature_vcf ? "--no-mature-suffix \"${stem}_nc_nomat_final.vcf\"" : ""
    """
    # Rewrite each segment FASTA header to exactly the sample_id so nextclade
    # seqName matches sample_id in the result JSON (enables per-sample demux).
    mkdir -p prepared
    for fasta_file in inputs/*.s.fasta; do
        sample_id=\$(basename "\$fasta_file" .s.fasta | sed -E 's/_[^_]+\$//')
        awk -v sample_id="\$sample_id" '/^>/{print ">"sample_id; next}{print}' "\$fasta_file" > "prepared/\${sample_id}.fasta"
    done
    cat prepared/*.fasta > batch.fasta

    nextclade run \\
        --input-dataset ${dataset_dir} \\
        --output-all nc_output \\
        --output-basename batch \\
        --include-reference \\
        --alignment-preset ${params.nextclade_alignment_preset} \\
        --jobs ${threads} \\
        batch.fasta

    nextclade_post_process.py \\
        --json    nc_output/batch.json \\
        --tsv     nc_output/batch.tsv \\
        --segment ${segment} \\
        --status-summary ${segment}_nextclade_status.tsv

    REFERENCE_NAME=\$(head -n1 ${dataset_dir}/reference.fasta | sed -E 's/>(\\S+).*/\\1/g')
    SUBTYPE=\$(head -n1 ${dataset_dir}/reference.fasta | sed -E 's/.+\\(\\S+\\((\\S+)\\)\\S*\\).+/\\1/g')
    if echo "\$SUBTYPE" | grep -q '^>'; then SUBTYPE="unknown"; fi

    nextclade_to_freebayes_vcf.py \\
        -i nc_output/batch.json \\
        -d ${dataset_dir} \\
        --suffix "${stem}_nc_final.vcf" \\
        ${no_mature_arg}
    """
}
