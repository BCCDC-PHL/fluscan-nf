// shared/select_db.nf - Reference bundle selection for FASTA samples without an explicit bundle

process select_db {

    tag { sample_id }

    conda "${projectDir}/environments/fluscan.yml"

    publishDir "${params.outdir}/${sample_id}/db_selection", pattern: "${sample_id}_db_selection*.tsv", mode: 'copy'

    input:
    tuple val(sample_id), path(query), path(blast_db)

    output:
    tuple val(sample_id), path("${sample_id}_db_selection.tsv"), emit: selection
    tuple val(sample_id), path("${sample_id}_db.txt"),           emit: db
    tuple val(sample_id), path("${sample_id}_db_selection_blastn.tsv")
    tuple val(sample_id), path("${sample_id}-select_db-provenance.yml"), emit: provenance

    script:
    blast_db_name = blast_db[0]
    """
    # Every subject is kept (no -max_target_seqs 1): select_db.py ranks the bundles itself.
    blastn \\
        -num_threads ${task.cpus} \\
        -max_hsps 1 \\
        -query ${query} \\
        -task blastn \\
        -db ${blast_db_name} \\
        -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore score" \\
        > ${sample_id}_db_selection_blastn.tsv

    select_db.py ${sample_id}_db_selection_blastn.tsv \\
        --query ${query} \\
        --sample-id ${sample_id} \\
        --min-identity ${params.db_select_min_identity} \\
        --min-coverage ${params.db_select_min_coverage} \\
        --out-tsv ${sample_id}_db_selection.tsv \\
        --out-db ${sample_id}_db.txt

    printf -- "- process_name: select_db\\n" > ${sample_id}-select_db-provenance.yml
    printf -- "  min_identity: ${params.db_select_min_identity}\\n  min_coverage: ${params.db_select_min_coverage}\\n" >> ${sample_id}-select_db-provenance.yml
    printf -- "  selected_db: '%s'\\n" "\$(cat ${sample_id}_db.txt)" >> ${sample_id}-select_db-provenance.yml
    printf -- "  tools: \\n  - tool_name: blastn\\n    tool_version: \$(blastn -version 2>&1 | head -n1 | cut -d' ' -f2)\\n" >> ${sample_id}-select_db-provenance.yml
    """
}
