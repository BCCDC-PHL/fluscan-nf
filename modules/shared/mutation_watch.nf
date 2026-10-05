// shared/mutation_watch.nf - Scan a sample's final annotated VCFs for custom watchlist substitutions.
//
// Runs in all three workflows off whatever each one calls its final VCFs; the only thing it
// needs is a VCF carrying ANN=. One task per sample scans every VCF against every watchlist.
//
// Publishes one *_watchlist_hits.tsv per (VCF, watchlist), with one row per matching variant.
// A scan that finds nothing writes a header-only file.
//
// Nothing here writes to the VCF. FLUMUT= is curated and always-on and belongs in the VCF;
// a watchlist is set on our own terms, and merging it in would make the same consensus
// produce different VCFs depending on the contents of a CSV.

// One scan per (VCF, watchlist): [vcf, watchlist, segment, output stem]. segments[i] labels
// vcfs[i]; the watchlist name is in the stem so several lists never collide on publish.
def watch_scans(sample_id, segments, vcfs, watchlists) {
    [segments, [vcfs].flatten()].transpose().collectMany { segment, vcf ->
        [watchlists].flatten().collect { watchlist ->
            def list_name = watchlist.baseName.replaceAll(/_watchlist$/, '')
            [vcf, watchlist, segment, "${sample_id}_${segment}_${list_name}"]
        }
    }
}

process mutation_watch {

    tag { sample_id }

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
    tuple val(sample_id), val(segments), path(vcfs, stageAs: 'vcf/*')
    path(watchlists, stageAs: 'watchlists/*')

    output:
    tuple val(sample_id), path("*_watchlist_hits.tsv"),                      emit: hits
    tuple val(sample_id), path("${sample_id}_mutation_watch_provenance.yml"), emit: provenance

    script:
    def scans = watch_scans(sample_id, segments, vcfs, watchlists).collect { vcf, watchlist, segment, stem ->
        "mutation_watcher.py --vcf ${vcf} --watchlist ${watchlist} --sample ${sample_id} --segment '${segment}' -o ${stem}_watchlist_hits.tsv"
    }
    """
    for watchlist in watchlists/*; do
        printf -- "- process_name: mutation_watch\\n"                                             >> ${sample_id}_mutation_watch_provenance.yml
        printf -- "  tools:\\n"                                                                   >> ${sample_id}_mutation_watch_provenance.yml
        printf -- "    - tool_name: python\\n      tool_version: \$(python3 --version 2>&1 | cut -d' ' -f2)\\n" >> ${sample_id}_mutation_watch_provenance.yml
        printf -- "    - tool_name: mutation_watcher.py\\n      tool_version: \$(mutation_watcher.py --version | awk '{print \$NF}')\\n" >> ${sample_id}_mutation_watch_provenance.yml
        printf -- "  watchlist: \$(basename "\$watchlist")\\n"                                    >> ${sample_id}_mutation_watch_provenance.yml
    done

    ${scans.join('\n    ')}
    """

    stub:
    def hits = watch_scans(sample_id, segments, vcfs, watchlists).collect { it[3] + '_watchlist_hits.tsv' }
    """
    touch ${hits.join(' ')}
    touch ${sample_id}_mutation_watch_provenance.yml
    """
}
