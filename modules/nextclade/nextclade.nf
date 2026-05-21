// nextclade/nextclade.nf - Per-segment Nextclade alignment and VCF conversion processes

process nextclade_run {

    conda "${projectDir}/environments/nextclade.yml"

    tag { "${sample_id}_${segment}" }

    publishDir "${params.outdir}/${sample_id}/${workflow}/json", pattern: "${sample_id}_${segment}_nextclade.json", mode: 'copy'
    publishDir "${params.outdir}/${sample_id}/${workflow}/tsv", pattern: "${sample_id}_${segment}_nextclade.tsv",  mode: 'copy'

    input:
    tuple val(sample_id), val(segment), path(consensus_fa), path(dataset_dir)

    output:
    tuple val(sample_id), val(segment), path("${sample_id}_${segment}_nextclade.json"), path(dataset_dir), emit: json_with_ref
    tuple val(sample_id), val(segment), path("${sample_id}_${segment}_nextclade.tsv"),                      emit: tsv

    script:
    workflow = task.ext.workflow ?: "unknown_workflow"

    """
    nextclade run \
        --input-dataset ${dataset_dir} \
        --output-json ${sample_id}_${segment}_nextclade.json \
        --output-tsv  ${sample_id}_${segment}_nextclade.tsv \
        ${consensus_fa}
    """
}


process nextclade_to_vcf {

    conda "${projectDir}/environments/nextclade.yml"

    tag { "${sample_id}_${segment}" }

    //publishDir "${params.outdir}/${sample_id}/${workflow}", pattern: "${sample_id}*nc.vcf", mode: 'copy'

    input:
    tuple val(sample_id), val(segment), path(nextclade_json), path(dataset_dir)

    output:
    tuple val(sample_id), val(segment), path("${sample_id}*nc.vcf"), path(dataset_dir), emit: vcf

    script:
    workflow = task.ext.workflow ?: "unknown_workflow"
    """
    REFERENCE_NAME=\$(head -n1 ${dataset_dir}/reference.fasta | sed -E 's/>(\\S+).+/\\1/g')
    SUBTYPE=\$(head -n1 ${dataset_dir}/reference.fasta | sed -E 's/.+\\(\\S+\\((\\S+)\\)\\S*\\).+/\\1/g')
    SEGMENT=${dataset_dir}
    SEGMENT=\${SEGMENT^^} # uppercase

    nextclade_to_freebayes_vcf.py \
        -i ${nextclade_json} \
        -d ${dataset_dir} \
        -o ${sample_id}__\${REFERENCE_NAME}__\${SUBTYPE}__\${SEGMENT}_nc.vcf \
        --sample ${sample_id}
    """
}


// Cross-sample per-segment batch: runs nextclade once for all samples' sequences
// for a given segment, then splits output into per-sample JSON and TSV files.
process nextclade_run_batch {

    conda "${projectDir}/environments/nextclade.yml"

    tag { segment }

    publishDir "${params.outdir}", mode: 'copy',
        saveAs: { filename ->
            if (filename =~ /.*_nextclade\.(json|tsv)$/) {
                def file_extension = filename.endsWith('.json') ? 'json' : 'tsv'
                def filename_base  = filename.replaceAll(/_nextclade\.(json|tsv)$/, '')
                def sample_id      = filename_base.substring(0, filename_base.lastIndexOf('_'))
                def workflow       = task.ext.workflow ?: 'unknown_workflow'
                return "${sample_id}/${workflow}/${file_extension}/${filename}"
            }
            return null
        }

    input:
    tuple val(segment), val(sample_ids),
          path(consensus_fastas, stageAs: 'inputs/*'),
          path(dataset_dir)

    output:
    tuple val(segment), val(sample_ids), path("*_nextclade.json"), path(dataset_dir), emit: json_with_ref
    tuple val(segment), path("*_nextclade.tsv"), emit: tsv

    script:
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
        --output-json   batch_nextclade.json \\
        --output-tsv    batch_nextclade.tsv \\
        batch.fasta

    # Split batch JSON into one file per sample
    python3 -c "
import json
batch_data = json.load(open('batch_nextclade.json'))
for result in batch_data.get('results', []):
    sample_id = result['seqName']
    sample_data = dict(batch_data)
    sample_data['results'] = [result]
    with open(sample_id + '_${segment}_nextclade.json', 'w') as json_out:
        json.dump(sample_data, json_out)
"

    # Split batch TSV into one file per sample (header + one data row each)
    python3 -c "
import csv
with open('batch_nextclade.tsv') as tsv_in:
    reader = csv.DictReader(tsv_in, delimiter='\\t')
    for row in reader:
        sample_id = row['seqName']
        with open(sample_id + '_${segment}_nextclade.tsv', 'w', newline='') as tsv_out:
            writer = csv.DictWriter(tsv_out, fieldnames=reader.fieldnames, delimiter='\\t')
            writer.writeheader()
            writer.writerow(row)
"
    """
}


// Converts per-sample nextclade JSON files (produced by nextclade_run_batch) into
// freebayes-style VCFs, one per sample, for the given segment.
process nextclade_to_vcf_batch {

    conda "${projectDir}/environments/nextclade.yml"

    tag { segment }

    input:
    tuple val(segment), val(sample_ids),
          path(json_files, stageAs: 'json/*'),
          path(dataset_dir)

    output:
    tuple val(segment), path('*_nc.vcf'), path(dataset_dir), emit: vcf

    script:
    """
    REFERENCE_NAME=\$(head -n1 ${dataset_dir}/reference.fasta | sed -E 's/>(\\S+).+/\\1/g')
    SUBTYPE=\$(head -n1 ${dataset_dir}/reference.fasta | sed -E 's/.+\\(\\S+\\((\\S+)\\)\\S*\\).+/\\1/g')
    SEGMENT="${segment}"

    for json_file in json/*_nextclade.json; do
        json_basename=\$(basename "\$json_file" _nextclade.json)
        sample_id=\${json_basename%_${segment}}
        nextclade_to_freebayes_vcf.py \\
            -i "\$json_file" \\
            -d ${dataset_dir} \\
            -o "\${sample_id}__\${REFERENCE_NAME}__\${SUBTYPE}__\${SEGMENT}_nc.vcf" \\
            --sample "\$sample_id"
    done
    """
}
