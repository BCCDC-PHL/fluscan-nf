// shared/mutation_watch.nf - Scan a final annotated VCF for a custom watchlist's substitutions.
//
// Runs in all three workflows off whatever each one calls its final VCF; the only thing it
// needs is a VCF carrying ANN=.
//
// Publishes one file, *_watchlist_hits.tsv, with one row per matching variant. A scan that
// finds nothing writes a header-only file.
//
// Nothing here writes to the VCF. FLUMUT= is curated and always-on and belongs in the VCF;
// a watchlist is set on our own terms, and merging it in would make the same consensus
// produce different VCFs depending on the contents of a CSV.

process mutation_watch {

    tag { segment ? "${sample_id}_${segment}_${watchlist.baseName}" : "${sample_id}_${watchlist.baseName}" }

    conda "${projectDir}/environments/fluscan.yml"

    publishDir "${params.outdir}", mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith('.tsv')) {
                def wf = task.ext.workflow ?: 'watchlist'
                return "${sample_id}/${wf}/${filename}"
            }
            // Provenance is not published from here: collect_provenance concatenates it
            // with the rest of the sample's records and publishes that single file.
            return null
        }

    input:
    tuple val(sample_id), val(segment), path(vcf), path(watchlist)

    output:
    tuple val(sample_id), val(segment), path("*_watchlist_hits.tsv"), emit: hits
    // Globbed rather than named: the stem is computed in the script block from the
    // watchlist's own filename, and an output declaration must not depend on that.
    tuple val(sample_id), path("*_mutation_watch_provenance.yml"), emit: provenance

    script:
    // The watchlist name is part of the output stem so several lists can be scanned for in
    // one run without colliding in the publish directory.
    def list_name = watchlist.baseName.replaceAll(/_watchlist$/, '')
    def prefix = segment ? "${segment}_${list_name}" : list_name
    def stem = "${sample_id}_${prefix}"
    """
    printf -- "- process_name: mutation_watch\\n"                                                 > ${stem}_mutation_watch_provenance.yml
    printf -- "  tools:\\n"                                                                       >> ${stem}_mutation_watch_provenance.yml
    printf -- "    - tool_name: python\\n      tool_version: \$(python3 --version 2>&1 | cut -d' ' -f2)\\n" >> ${stem}_mutation_watch_provenance.yml
    printf -- "    - tool_name: mutation_watcher.py\\n      tool_version: \$(mutation_watcher.py --version | awk '{print \$NF}')\\n" >> ${stem}_mutation_watch_provenance.yml
    printf -- "  watchlist: ${watchlist.name}\\n"                                                 >> ${stem}_mutation_watch_provenance.yml

    mutation_watcher.py \\
        --vcf ${vcf} \\
        --watchlist ${watchlist} \\
        --sample ${sample_id} \\
        --segment '${segment}' \\
        -o ${stem}_watchlist_hits.tsv
    """

    stub:
    def list_name = watchlist.baseName.replaceAll(/_watchlist$/, '')
    def prefix = segment ? "${segment}_${list_name}" : list_name
    def stem = "${sample_id}_${prefix}"
    """
    touch ${stem}_watchlist_hits.tsv
    touch ${stem}_mutation_watch_provenance.yml
    """
}
