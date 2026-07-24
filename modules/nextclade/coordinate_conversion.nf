// nextclade/coordinate_conversion.nf - Build a source->target reference coordinate map.
//
// Runs once per pipeline invocation (not per-sample). Translates nucleotide and
// amino-acid indices from flumut's fixed reference (the source) onto the selected
// Nextclade reference dataset (the target) coordinate system, so a downstream step
// can remap flumut mutations into Nextclade coordinates. Idempotent: given the same
// two reference directories it always produces the same map, including the trivial
// identity case where the Nextclade dataset is derived from flumut's reference.
//
// This step is independent of, and does not modify, nextclade_to_freebayes_vcf.py
// or any of its consuming processes (nextclade_to_vcf*, post_process_vcf*).

process coordinate_conversion {

    tag { target_name }

    conda "${projectDir}/environments/fluscan.yml"

    publishDir "${params.outdir}/aggregate", pattern: "reference_remap", mode: 'copy'

    input:
    tuple val(source_name), path(source_dataset_dir, stageAs: 'source_dataset'),
          val(target_name), path(target_dataset_dir, stageAs: 'target_dataset')

    output:
    path("reference_remap"),                                          emit: map_dir
    path("reference_remap/coordinate_map.json"),                      emit: json
    path("reference_remap/coordinate_map.*.json"),                    emit: json_shards
    path("reference_remap/coordinate_conversion_provenance.yml"),     emit: provenance

    script:
    """
    mkdir -p reference_remap

    printf -- "- process_name: coordinate_conversion\\n"                                              > reference_remap/coordinate_conversion_provenance.yml
    printf -- "  tools:\\n"                                                                          >> reference_remap/coordinate_conversion_provenance.yml
    printf -- "    - tool_name: python\\n      tool_version: \$(python3 --version 2>&1 | cut -d' ' -f2)\\n" >> reference_remap/coordinate_conversion_provenance.yml
    printf -- "    - tool_name: biopython\\n      tool_version: \$(python3 -c 'import Bio; print(Bio.__version__)')\\n" >> reference_remap/coordinate_conversion_provenance.yml
    printf -- "    - tool_name: build_reference_coordinate_map.py\\n      tool_version: 2.2.0\\n"     >> reference_remap/coordinate_conversion_provenance.yml
    printf -- "  parameters:\\n"                                                                     >> reference_remap/coordinate_conversion_provenance.yml
    printf -- "    - parameter: source_reference\\n      value: ${source_name}\\n"                   >> reference_remap/coordinate_conversion_provenance.yml
    printf -- "    - parameter: target_reference\\n      value: ${target_name}\\n"                   >> reference_remap/coordinate_conversion_provenance.yml

    build_reference_coordinate_map.py \\
        --source-ref-dir source_dataset \\
        --target-ref-dir target_dataset \\
        --out-dir reference_remap \\
        --source-name ${source_name} \\
        --target-name ${target_name}
    """

    stub:
    """
    mkdir -p reference_remap
    touch reference_remap/coordinate_map.json
    for seg in PB2 PB1 PA HA NP NA M NS; do
        touch reference_remap/coordinate_map.\${seg}.json
    done
    touch reference_remap/coordinate_map.nt.tsv
    touch reference_remap/coordinate_map.aa.tsv
    touch reference_remap/coordinate_map.qc.tsv
    touch reference_remap/${source_name}_alignment.fasta
    touch reference_remap/${source_name}_aa_alignment.fasta
    touch reference_remap/coordinate_conversion_provenance.yml
    """
}
