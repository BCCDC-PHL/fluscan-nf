// nextclade/merge_flumut_mutations.nf - Annotate the Nextclade VCF with FluMutDB markers.
//
// Consumes the flumut->Nextclade coordinate map (per-segment shard) and a per-sample/
// per-segment Nextclade freebayes-style VCF. For every AA change already annotated in a
// row's ANN= field, merge_flumut_mutations.py queries FluMutDB (shipped by the flumutdb
// conda package) for a matching curated marker and, on a hit, adds a FLUMUT= INFO field.
// No rows are added, removed or reordered. Writes a NEW *_flumut_merged.vcf; the input
// *_final.vcf is untouched.
//
// The marker set comes from the database, so no per-sample flumut output is needed here
// (the flumut process still runs and publishes its own TSVs).
//
// This step does not modify any existing process or script (nextclade_to_freebayes_vcf.py,
// post_process_vcf*, flumut, coordinate_conversion).

process merge_flumut_mutations {

    tag { "${sample_id}_${segment}" }

    conda "${projectDir}/environments/fluscan.yml"

    publishDir "${params.outdir}", mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith('_flumut_merged.vcf')) {
                def wf = task.ext.workflow ?: 'nextclade'
                return "${sample_id}/${wf}/${filename}"
            }
            // Provenance is not published from here: collect_provenance concatenates it
            // with the rest of the sample's records and publishes that single file.
            return null
        }

    input:
    tuple val(sample_id), val(segment), path(vcf), path(coord_json)

    output:
    tuple val(sample_id), val(segment), path("*_flumut_merged.vcf"),        emit: vcf
    tuple val(sample_id), path("${sample_id}_${segment}_merge_flumut_provenance.yml"), emit: provenance

    script:
    // Preserve the sample__ref__subtype__segment stem from the input VCF when present.
    def stem = vcf.name.replaceAll(/(_final)?\.vcf$/, '')
    """
    printf -- "- process_name: merge_flumut_mutations\\n"                                         > ${sample_id}_${segment}_merge_flumut_provenance.yml
    printf -- "  tools:\\n"                                                                       >> ${sample_id}_${segment}_merge_flumut_provenance.yml
    printf -- "    - tool_name: python\\n      tool_version: \$(python3 --version 2>&1 | cut -d' ' -f2)\\n" >> ${sample_id}_${segment}_merge_flumut_provenance.yml
    printf -- "    - tool_name: merge_flumut_mutations.py\\n      tool_version: \$(merge_flumut_mutations.py --version | awk '{print \$NF}')\\n" >> ${sample_id}_${segment}_merge_flumut_provenance.yml
    printf -- "    - tool_name: flumutdb\\n      tool_version: \$(python3 -c 'from flumut.DbReader import get_db_version; print("%s.%s" % get_db_version()[:2])')\\n" >> ${sample_id}_${segment}_merge_flumut_provenance.yml

    merge_flumut_mutations.py \\
        --vcf ${vcf} \\
        --coord-map ${coord_json} \\
        --segment ${segment} \\
        --sample ${sample_id} \\
        -o ${stem}_flumut_merged.vcf
    """

    stub:
    def stem = vcf.name.replaceAll(/(_final)?\.vcf$/, '')
    """
    touch ${stem}_flumut_merged.vcf
    touch ${sample_id}_${segment}_merge_flumut_provenance.yml
    """
}
