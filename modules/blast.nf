// blast.nf - BLAST-related processes

process makeblastdb {

    tag { fasta_db.name }

    input:
    path(fasta_db)

    output:
    tuple path("${fasta_db.name}"), path("${fasta_db.name}.*"), emit: db

    script:
    """
    makeblastdb -in ${fasta_db} -dbtype nucl
    """
}


process blastn {

    errorStrategy 'ignore'

    label 'medium'

    tag {sample_id}

    publishDir "${params.outdir}/${sample_id}/", pattern: "${sample_id}*blastn.tsv" , mode:'copy'
    publishDir "${params.outdir}/${sample_id}/", pattern: "${sample_id}*filter.tsv" , mode:'copy'
    publishDir "${params.outdir}/${sample_id}/", pattern: "${sample_id}*fasta" , mode:'copy'

    input: 
    tuple val(sample_id), path(query)
    tuple path(blast_db_fasta), path(blast_db_index_files)

    output:
    tuple val(sample_id), path("${sample_id}*blastn.tsv"), emit: raw
    tuple val(sample_id), path("${sample_id}*filter.tsv"), path("${sample_id}*fasta"), emit: main
    tuple val(sample_id), path("${sample_id}*_ref.fasta"), emit: ref
    path("${sample_id}*filter.tsv"), emit: filter
    tuple val(sample_id), path("${sample_id}-*-provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: blastn\\n" > ${sample_id}-blastn-provenance.yml
    printf -- "  tools: \\n  - tool_name: blastn\\n    tool_version: \$(blastn -version 2>&1 | head -n1 | cut -d' ' -f2)\\n" >> ${sample_id}-blastn-provenance.yml

    blastn -num_threads ${task.cpus} -query ${query} -db ${blast_db_fasta} -outfmt 6 > ${sample_id}_blastn.tsv &&

    filter_alignments.py ${sample_id}_blastn.tsv \
    --seqs ${blast_db_fasta} \
    --out_tsv ${sample_id}_blastn_filter.tsv \
    --out_fasta ${sample_id}_ref.fasta 
    """
}
