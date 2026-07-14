// nextclade/insert_flumut_mutations.nf - Merge flumut AA markers into the Nextclade VCF.
//
// Consumes the (already-built but previously orphaned) flumut->Nextclade coordinate map
// (per-segment shard), flumut's per-sample mutations TSV, and a per-sample/per-segment
// Nextclade freebayes-style VCF. Labels every AA change with DETECTED_BY=fluscan|flumut|both
// and appends new rows for flumut-only markers. Writes a NEW *_flumut_merged.vcf; the input
// *_final.vcf is untouched. Target CDS placement comes from the map's target_cds block, so
// no GFF/dataset dir is needed here.
//
// This step does not modify any existing process or script (nextclade_to_freebayes_vcf.py,
// post_process_vcf*, flumut, coordinate_conversion).

process insert_flumut_mutations {

    tag { "${sample_id}_${segment}" }

    conda "${projectDir}/environments/fluscan.yml"

    publishDir "${params.outdir}", mode: 'copy',
        saveAs: { filename ->
            if (filename.endsWith('_flumut_merged.vcf')) {
                def wf = task.ext.workflow ?: 'nextclade'
                return "${sample_id}/${wf}/${filename}"
            }
            if (filename.endsWith('_provenance.yml')) {
                return null
            }
            return null
        }

    input:
    tuple val(sample_id), val(segment), path(vcf), path(coord_json),
          path(flumut_mutations), path(consensus)

    output:
    tuple val(sample_id), val(segment), path("*_flumut_merged.vcf"),        emit: vcf
    tuple val(sample_id), path("${sample_id}_${segment}_insert_flumut_provenance.yml"), emit: provenance

    script:
    // Preserve the sample__ref__subtype__segment stem from the input VCF when present.
    def stem = vcf.name.replaceAll(/(_final)?\.vcf$/, '')
    """
    printf -- "- process_name: insert_flumut_mutations\\n"                                        > ${sample_id}_${segment}_insert_flumut_provenance.yml
    printf -- "  tools:\\n"                                                                       >> ${sample_id}_${segment}_insert_flumut_provenance.yml
    printf -- "    - tool_name: python\\n      tool_version: \$(python3 --version 2>&1 | cut -d' ' -f2)\\n" >> ${sample_id}_${segment}_insert_flumut_provenance.yml
    printf -- "    - tool_name: insert_flumut_mutations.py\\n      tool_version: 1.0.0\\n"         >> ${sample_id}_${segment}_insert_flumut_provenance.yml

    insert_flumut_mutations.py \\
        --vcf ${vcf} \\
        --coord-map ${coord_json} \\
        --flumut-mutations ${flumut_mutations} \\
        --consensus ${consensus} \\
        --segment ${segment} \\
        --sample ${sample_id} \\
        -o ${stem}_flumut_merged.vcf
    """

    stub:
    def stem = vcf.name.replaceAll(/(_final)?\.vcf$/, '')
    """
    touch ${stem}_flumut_merged.vcf
    touch ${sample_id}_${segment}_insert_flumut_provenance.yml
    """
}
