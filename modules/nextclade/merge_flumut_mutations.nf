// nextclade/merge_flumut_mutations.nf - Annotate the Nextclade VCFs with FluMutDB markers.
//
// Consumes a per-segment coordinate-map shard and every sample's VCF for that segment, and
// adds a FLUMUT= INFO field to rows whose ANN= AA change matches a curated FluMutDB marker. No
// rows are added, removed or reordered. Writes a NEW *_flumut_annotated.vcf per input, leaving
// each *_final.vcf untouched. Markers come from the flumutdb package, not a per-sample run.

process merge_flumut_mutations {

    tag { "${db_name}/${segment}" }

    conda "${projectDir}/environments/fluscan.yml"

    errorStrategy { task.exitStatus in 1..140 ? 'retry' : 'terminate' }
    maxRetries 2

    publishDir "${params.outdir}", mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith('_flumut_annotated.vcf')) {
                def sample_id = filename.split('__')[0]
                return "${sample_id}/nextclade/flumut_annotated/${filename}"
            }
            // Provenance is not published from here: collect_provenance concatenates it
            // with the rest of each sample's records and publishes that single file.
            return null
        }

    input:
    tuple val(db_name), val(segment), path(vcfs, stageAs: 'vcf/*'), path(coord_json)

    output:
    tuple val(db_name), val(segment), path("*_flumut_annotated.vcf"),                  emit: vcf
    tuple val(db_name), val(segment), path("${segment}_merge_flumut_provenance.yml"),  emit: provenance

    script:
    """
    printf -- "- process_name: merge_flumut_mutations\\n"                                         > ${segment}_merge_flumut_provenance.yml
    printf -- "  tools:\\n"                                                                       >> ${segment}_merge_flumut_provenance.yml
    printf -- "    - tool_name: python\\n      tool_version: \$(python3 --version 2>&1 | cut -d' ' -f2)\\n" >> ${segment}_merge_flumut_provenance.yml
    printf -- "    - tool_name: merge_flumut_mutations.py\\n      tool_version: \$(merge_flumut_mutations.py --version | awk '{print \$NF}')\\n" >> ${segment}_merge_flumut_provenance.yml
    printf -- "    - tool_name: flumutdb\\n      tool_version: \$(python3 -c 'from flumut.DbReader import get_db_version; print("%s.%s" % get_db_version()[:2])')\\n" >> ${segment}_merge_flumut_provenance.yml

    merge_flumut_mutations.py --vcf vcf/*.vcf --coord-map ${coord_json} --segment ${segment}
    """

    stub:
    """
    for vcf in vcf/*.vcf; do
        stem=\$(basename "\$vcf" .vcf)
        touch "\${stem%_final}_flumut_annotated.vcf"
    done
    touch ${segment}_merge_flumut_provenance.yml
    """
}
