// shared/mutation_watch.nf - Scan a sample's final annotated VCFs for custom watchlist substitutions.
//
// Runs in all three workflows off whatever each one calls its final VCFs; the only thing it
// needs is a VCF carrying ANN=. One task per sample scans every VCF against every watchlist.
//
// Publishes one watchlist/<list>/*_watchlist_hits.tsv per (VCF, watchlist), with one row per
// matching variant. A scan that finds nothing writes a header-only file.
//
// Nothing here writes to the VCF. FLUMUT= is curated and always-on and belongs in the VCF;
// a watchlist is set on our own terms, and merging it in would make the same consensus
// produce different VCFs depending on the contents of a CSV.

// One scan per (VCF, watchlist): [vcf, watchlist, segment, output path]. segments[i] labels
// vcfs[i]; the workflow tag keeps FASTA mode's PAFTOOLS and NEXTCLADE hits apart in one list dir.
def watch_scans(sample_id, workflow, segments, vcfs, watchlists) {
    def wf_tag = [paftools: 'pf', nextclade: 'nc', short_read: 'sr'][workflow]
    [segments, [vcfs].flatten()].transpose().collectMany { segment, vcf ->
        [watchlists].flatten().collect { watchlist ->
            def list_name = watchlist.baseName.replaceAll(/_watchlist$/, '')
            [vcf, watchlist, segment, "${list_name}/${sample_id}_${segment}_${list_name}_${wf_tag}_watchlist_hits.tsv"]
        }
    }
}

process mutation_watch {

    tag { sample_id }

    conda "${projectDir}/environments/fluscan.yml"

    // Provenance is not published from here: collect_provenance concatenates it with the
    // rest of the sample's records and publishes that single file.
    publishDir "${params.outdir}/${sample_id}/watchlist", pattern: "*/*_watchlist_hits.tsv", mode: 'copy'

    input:
    tuple val(sample_id), val(segments), path(vcfs, stageAs: 'vcf/*')
    path(watchlists, stageAs: 'watchlists/*')

    output:
    tuple val(sample_id), path("*/*_watchlist_hits.tsv"),                    emit: hits
    tuple val(sample_id), path("${sample_id}_mutation_watch_provenance.yml"), emit: provenance

    script:
    def scans = watch_scans(sample_id, task.ext.workflow, segments, vcfs, watchlists)
    def commands = scans.collect { vcf, watchlist, segment, out ->
        "mutation_watcher.py --vcf ${vcf} --watchlist ${watchlist} --sample ${sample_id} --segment '${segment}' -o ${out}"
    }
    """
    for watchlist in watchlists/*; do
        printf -- "- process_name: mutation_watch\\n"                                             >> ${sample_id}_mutation_watch_provenance.yml
        printf -- "  tools:\\n"                                                                   >> ${sample_id}_mutation_watch_provenance.yml
        printf -- "    - tool_name: python\\n      tool_version: \$(python3 --version 2>&1 | cut -d' ' -f2)\\n" >> ${sample_id}_mutation_watch_provenance.yml
        printf -- "    - tool_name: mutation_watcher.py\\n      tool_version: \$(mutation_watcher.py --version | awk '{print \$NF}')\\n" >> ${sample_id}_mutation_watch_provenance.yml
        printf -- "  watchlist: \$(basename "\$watchlist")\\n"                                    >> ${sample_id}_mutation_watch_provenance.yml
    done

    mkdir -p ${scans.collect { it[3].tokenize('/')[0] }.unique().join(' ')}
    ${commands.join('\n    ')}
    """

    stub:
    def scans = watch_scans(sample_id, task.ext.workflow, segments, vcfs, watchlists)
    """
    mkdir -p ${scans.collect { it[3].tokenize('/')[0] }.unique().join(' ')}
    touch ${scans.collect { it[3] }.join(' ')}
    touch ${sample_id}_mutation_watch_provenance.yml
    """
}
