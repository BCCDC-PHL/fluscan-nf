// paftools/paftools.nf - per-sample consensus-to-reference variant calling for the PAFTOOLS workflow

// Every segment of one sample in one task: minimap2 alignment, coverage QC, paftools.js
// call, SnpEff, freebayes-style conversion and MNV codon correction. A segment with no
// BLAST reference, or that minimap2 cannot align, gets a FAIL coverage row and no VCF.
process paftools_variants {

    tag { sample_id }

    conda "${projectDir}/environments/fluscan.yml"

    publishDir "${params.outdir}/${sample_id}/qc",       pattern: "*_coverage_qc.tsv", mode: 'copy'
    publishDir "${params.outdir}/${sample_id}/paftools", pattern: "*_final.vcf",       mode: 'copy'

    input:
    tuple val(sample_id), path(consensus), path(reference), path(snpeff_dir, stageAs: 'snpeff')

    output:
    tuple val(sample_id), path("${sample_id}__*_paf_final.vcf"),                emit: vcf
    tuple val(sample_id), path("${sample_id}_*_coverage_qc.tsv"),               emit: coverage
    tuple val(sample_id), path("${sample_id}_paftools_variants_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: paftools_variants\\n" > ${sample_id}_paftools_variants_provenance.yml
    printf -- "  tools: \\n" >> ${sample_id}_paftools_variants_provenance.yml
    printf -- "  - tool_name: minimap2\\n    tool_version: \$(minimap2 --version 2>&1 | head -n1)\\n" >> ${sample_id}_paftools_variants_provenance.yml
    printf -- "  - tool_name: paftools.js\\n    tool_version: \$(paftools.js version)\\n" >> ${sample_id}_paftools_variants_provenance.yml
    printf -- "  - tool_name: snpeff\\n    tool_version: \$(snpEff -version 2>&1 | cut -d\$'\\t' -f2)\\n" >> ${sample_id}_paftools_variants_provenance.yml
    printf -- "  - tool_name: python\\n    tool_version: \$(python --version 2>&1 | cut -d' ' -f2)\\n" >> ${sample_id}_paftools_variants_provenance.yml
    printf -- "  - tool_name: fix_mnv_codons.py\\n    tool_version: 1.0.0\\n" >> ${sample_id}_paftools_variants_provenance.yml

    # Reference headers are accession|subtype|segment; each per-segment copy keeps only the
    # accession, which is the chromosome name in the SnpEff database.
    awk -F'|' '/^>/{out = "ref_" \$NF ".fa"; print \$1 > out; next} {print > out}' ${reference}
    split_fasta.py split ${consensus}

    coverage_qc() {
        consensus_coverage_qc.py \\
          --paf "\$1" \\
          --reference "\$2" \\
          --sample-id ${sample_id} \\
          --segment "\$segment" \\
          --pass-pct ${params.consensus_coverage_pass_pct} \\
          --warn-pct ${params.consensus_coverage_warn_pct} \\
          --out-tsv ${sample_id}_\${segment}_coverage_qc.tsv
    }

    for consensus_fa in ${sample_id}_*.s.fasta; do
        segment=\$(basename "\$consensus_fa" .s.fasta)
        segment="\${segment#${sample_id}_}"
        reference_fa="ref_\${segment}.fa"

        [ -s "\$reference_fa" ] || { echo "WARN: no reference for \$segment" >&2; coverage_qc /dev/null /dev/null; continue; }
        minimap2 -cx asm20 --cs "\$reference_fa" "\$consensus_fa" > "\${segment}_raw.paf" \\
            || { echo "WARN: minimap2 failed for \$segment" >&2; coverage_qc /dev/null "\$reference_fa"; continue; }
        sort -k6,6 -k8,8n "\${segment}_raw.paf" > "\${segment}.paf"
        coverage_qc "\${segment}.paf" "\$reference_fa"

        IFS='|' read -r accession subtype _ <<< "\$(grep -m1 "|\${segment}\$" ${reference} | cut -c2-)"
        stem="${sample_id}__\${accession}__\${subtype}__\${segment}"

        paftools.js call -l 100 -L 100 -f "\$reference_fa" "\${segment}.paf" > "\${segment}_call.vcf"

        # -dataDir overrides the config's data.dir, which points wherever the bundle was built.
        snpEff -Xmx4g -noStats -hgvs1LetterAa \\
            -config snpeff/snpEff.config \\
            -dataDir "\$PWD/snpeff" \\
            "influenza_\${subtype,,}" "\${segment}_call.vcf" > "\${segment}_snpeff.vcf"

        paftools_to_freebayes_vcf.py \\
            --input "\${segment}_snpeff.vcf" \\
            --coverage-tsv ${sample_id}_\${segment}_coverage_qc.tsv \\
            --output "\${stem}_paf.vcf"

        fix_mnv_codons.py \\
            --input "\${stem}_paf.vcf" \\
            --reference "\$reference_fa" \\
            --output-vcf "\${stem}_paf_final.vcf" \\
            --min-af ${params.mnv_min_af}
    done
    """
}
