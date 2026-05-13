// blast.nf - BLAST-related processes

process makeblastdb {

    tag { fasta_db.name }

    conda "${projectDir}/environments/fluscan.yml"

    input:
    path(fasta_db)

    output:
    path("${blast_db_name}*"), emit: db

    script:
    blast_db_name = "${task.ext.workflow}_blastn_db.fasta"
    """
    makeblastdb -dbtype nucl -in ${fasta_db} -out ${blast_db_name}
    cp ${fasta_db} ${blast_db_name}
    """
}


process blastn {

    errorStrategy 'ignore'

    label 'medium'

    tag {sample_id}

    conda "${projectDir}/environments/fluscan.yml"

    publishDir "${params.outdir}/${sample_id}/", pattern: "${sample_id}*blastn.tsv" , mode:'copy'
    publishDir "${params.outdir}/${sample_id}/", pattern: "${sample_id}*filter.tsv" , mode:'copy'
    publishDir "${params.outdir}/${sample_id}/", pattern: "${sample_id}*fasta" , mode:'copy'

    input: 
    tuple val(sample_id), path(query)
    path(blast_db)

    output:
    tuple val(sample_id), path("${sample_id}*blastn.tsv"), emit: raw
    tuple val(sample_id), path("${sample_id}*filter.tsv"), path("${sample_id}*fasta"), emit: main
    tuple val(sample_id), path("${sample_id}*_ref.fasta"), emit: ref
    path("${sample_id}*filter.tsv"), emit: filter
    tuple val(sample_id), path("${sample_id}-*-provenance.yml"), emit: provenance

    script:
    blast_db_name = blast_db[0]
    """
    printf -- "- process_name: blastn\\n" > ${sample_id}-blastn-provenance.yml
    printf -- "  tools: \\n  - tool_name: blastn\\n    tool_version: \$(blastn -version 2>&1 | head -n1 | cut -d' ' -f2)\\n" >> ${sample_id}-blastn-provenance.yml

    blastn \
        -num_threads ${task.cpus} \
        -max_target_seqs 1 \
        -query ${query} \
        -perc_identity ${params.blastn_perc_identity} \
        -db ${blast_db_name} \
        -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore score" \
        > ${sample_id}_blastn.tsv &&

    filter_alignments.py ${sample_id}_blastn.tsv \
    --seqs ${blast_db_name} \
    --out_tsv ${sample_id}_blastn_filter.tsv \
    --out_fasta ${sample_id}_ref.fasta 
    """
}
