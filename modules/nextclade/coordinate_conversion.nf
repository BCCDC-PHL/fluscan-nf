// nextclade/coordinate_conversion.nf - Build a FluMut->Nextclade reference coordinate map.
//
// Runs once per bundle, not per sample. Translates nucleotide and amino-acid
// indices from flumut's fixed reference (source, read from the installed flumutdb package)
// onto the selected Nextclade dataset (target), so merge_flumut_mutations can place FluMutDB
// positions in Nextclade coordinates. Idempotent: the same flumutdb and dataset always give
// the same map.

process coordinate_conversion {

    tag { target_name }

    conda "${projectDir}/environments/fluscan.yml"

    publishDir "${params.outdir}/aggregate/reference_remap", pattern: "reference_remap", mode: 'copy',
        saveAs: { db_name }

    input:
    tuple val(db_name), val(target_name), path(target_dataset_dir, stageAs: 'target_dataset')

    output:
    tuple val(db_name), path("reference_remap"),                                          emit: map_dir
    tuple val(db_name), path("reference_remap/coordinate_map.json"),                      emit: json
    tuple val(db_name), path("reference_remap/coordinate_map.*.json"),                    emit: json_shards
    tuple val(db_name), path("reference_remap/coordinate_conversion_provenance.yml"),     emit: provenance

    script:
    """
    mkdir -p reference_remap

    printf -- "- process_name: coordinate_conversion\\n"                                              > reference_remap/coordinate_conversion_provenance.yml
    printf -- "  tools:\\n"                                                                          >> reference_remap/coordinate_conversion_provenance.yml
    printf -- "    - tool_name: python\\n      tool_version: \$(python3 --version 2>&1 | cut -d' ' -f2)\\n" >> reference_remap/coordinate_conversion_provenance.yml
    printf -- "    - tool_name: biopython\\n      tool_version: \$(python3 -c 'import Bio; print(Bio.__version__)')\\n" >> reference_remap/coordinate_conversion_provenance.yml
    printf -- "    - tool_name: build_reference_coordinate_map.py\\n      tool_version: \$(build_reference_coordinate_map.py --version | awk '{print \$NF}')\\n" >> reference_remap/coordinate_conversion_provenance.yml
    printf -- "    - tool_name: flumutdb\\n      tool_version: \$(python3 -c 'from flumut.DbReader import get_db_version; print("%s.%s" % get_db_version()[:2])')\\n" >> reference_remap/coordinate_conversion_provenance.yml
    printf -- "  parameters:\\n"                                                                     >> reference_remap/coordinate_conversion_provenance.yml
    printf -- "    - parameter: target_reference\\n      value: ${target_name}\\n"                   >> reference_remap/coordinate_conversion_provenance.yml

    build_reference_coordinate_map.py \\
        --target-ref-dir target_dataset \\
        --out-dir reference_remap \\
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
    touch reference_remap/flumut_alignment.fasta
    touch reference_remap/flumut_aa_alignment.fasta
    touch reference_remap/coordinate_conversion_provenance.yml
    """
}
