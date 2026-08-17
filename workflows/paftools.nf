#!/usr/bin/env nextflow

/*
 *   PAFTOOLS workflow (consensus FASTA input)
 *   -----------------------------------------
 *   Processes pre-assembled consensus FASTA files through quality control,
 *   subtype-aware reference selection via BLAST, minimap2/paftools variant
 *   calling, bcftools VCF sorting, SnpEff annotation, and summary extraction.
 *
 *   The sample name comes from the input filename (modules/shared/sample_id.nf); input
 *   headers supply only the segment label. Accepted header formats:
 *     Underscore-delimited (legacy): >SAMPLE_ID_segmentN_SEGMENT
 *       e.g. >WIN-AH-2024-FAV-0415-1-OS_segment1_PB2
 *     Pipe-delimited (GISAID-style): >EPI_ISL_ACCESSION|subtype|SEGMENT|...|clade_SEGMENT
 *       e.g. >EPI_ISL_20399898|A_/_H5N1|NA|6|EPI5261179|2.3.4.4b_NA
 *   Cleaned headers are normalised to {sample_id}_{SEGMENT} by prepare_consensus.
 *
 *   Reference database header format (<--db>/blastn_db.fa):
 *     >accession|subtype|segment
 *   e.g.
 *     >NC_026438.1|H1N1|PB2
 *
 */

nextflow.enable.dsl = 2

include { hash_files }                                   from '../modules/shared/hash_files.nf'
include { pipeline_provenance }                          from '../modules/shared/provenance.nf'
include { collect_provenance }                           from '../modules/shared/provenance.nf'
include { prepare_consensus }                            from '../modules/shared/consensus_qc.nf'
include { split_fasta as split_fasta_reference }         from '../modules/shared/consensus_qc.nf'
include { split_fasta as split_fasta_consensus }         from '../modules/shared/consensus_qc.nf'
include { paftools_align }                               from '../modules/paftools/alignment.nf'
include { consensus_coverage_qc }                        from '../modules/paftools/alignment.nf'
include { makeblastdb }                                  from '../modules/shared/blast.nf'
include { blastn }                                       from '../modules/shared/blast.nf'
include { flumut }                                       from '../modules/shared/flumut.nf'
include { paftools_call }                                from '../modules/paftools/variant_calling.nf'
include { paftools_sort_vcf }                            from '../modules/paftools/variant_calling.nf'
include { paftools_snpeff_annotation }                   from '../modules/paftools/variant_calling.nf'
include { paftools_to_freebayes_vcf }                    from '../modules/paftools/variant_calling.nf'
include { post_process_vcf }                             from '../modules/shared/post_process.nf'
include { mutation_watch }                               from '../modules/shared/mutation_watch.nf'
include { resolve_db }                                   from '../modules/shared/database.nf'
include { consensus_fasta_channel }                      from '../modules/shared/sample_id.nf'



workflow PAFTOOLS {

    // Reference selection is BLAST-driven and every called variant is SnpEff-annotated,
    // so both components are mandatory for this branch.
    def db = resolve_db(params.db, ['blastn', 'snpeff'])

    ch_workflow_metadata = Channel.value([
        workflow.sessionId,
        workflow.runName,
        workflow.manifest.name,
        workflow.manifest.version,
        workflow.start,
        [
            database_name:    db.name,
            database_version: db.version,
            database_path:    db.root,
        ],
    ])

    ch_pipeline_provenance = pipeline_provenance(ch_workflow_metadata)

    // Discover consensus FASTA files in the input directory.
    // Each file is expected to contain 8 segments for a single sample.
    ch_fasta_input = consensus_fasta_channel(params.fasta_input)

    ch_blastn_db = Channel.fromPath(db.blastn)

    main:

        // The whole snpeff/ directory is staged so the config and the database it names
        // travel together; the process overrides data.dir with -dataDir.
        ch_snpeff_dir = Channel.fromPath(db.snpeff_dir, checkIfExists: true)

        // Provenance channel starts with just the sample IDs; per-process YAMLs are joined
        // on as they are produced, and concatenated at the end by collect_provenance.
        ch_provenance = ch_fasta_input.map{ it -> it[0] }

        // Hash the input FASTA. One file per sample, so no list-wrapping here.
        hash_files(ch_fasta_input.combine(Channel.of("fasta_input")))

        // --- Step 1: Validate and clean input consensus FASTA ---
        prepare_consensus(ch_fasta_input)

        // One run-level table of the per-segment consensus QC rows. The name carries the
        // workflow because main.nf runs both FASTA branches into the same outdir, and two
        // collectFile operations must not target one path.
        prepare_consensus.out.qc
            .map{ it -> it[1] }
            .collectFile(
                name: 'paftools_consensus_qc_summary.tsv',
                storeDir: "${params.outdir}/aggregate",
                keepHeader: true,                   // keep header from first file only
                skip: 1,                            // drop the header line from the rest
                sort: false
            )

        // --- Step 2: Build BLAST database (once, shared across all samples) ---
        makeblastdb(ch_blastn_db).first().set{ch_blastdb}

        // --- Step 3: BLAST cleaned consensus against reference database ---
        blastn(prepare_consensus.out.cleaned, ch_blastdb)

        // --- Step 3b: HPAI mutation screening (per sample, parallel branch) ---
        flumut(prepare_consensus.out.cleaned)

        ch_consensus_segments = split_fasta_consensus(prepare_consensus.out.cleaned)
        ch_reference_segments = split_fasta_reference(blastn.out.ref)

        ch_segment_sequences = ch_consensus_segments.join(ch_reference_segments)

        paftools_align(ch_segment_sequences)

        // --- Step 4b: Check reference coverage of each alignment ---
        consensus_coverage_qc(paftools_align.out.paf)

        consensus_coverage_qc.out.coverage
            .map{ it -> it[1] }
            .flatten()                              // one Path per emission -> avoids the [name, content...] branch
            .collectFile(
                name: 'paftools_qc_summary.tsv',
                storeDir: "${params.outdir}/aggregate",
                keepHeader: true,                   // keep header from first file only
                skip: 1,                            // drop the header line from the rest
                sort: false
            )

        // --- Step 5: Align consensus to reference and call variants with paftools ---
        paftools_call(paftools_align.out.paf)

        // --- Step 6: Sort VCF for downstream annotation ---
        paftools_sort_vcf(paftools_call.out.vcf)

        // --- Step 7: Annotate amino-acid impacts with SnpEff ---
        paftools_snpeff_annotation(paftools_sort_vcf.out.vcf.combine(ch_snpeff_dir))

        paftools_to_freebayes_vcf(
            paftools_snpeff_annotation.out.vcf.join(consensus_coverage_qc.out.coverage)
        )

        // --- Step 8: Post-process annotated VCF (MNV codon correction) ---
        post_process_vcf(
            paftools_to_freebayes_vcf.out.vcf.join(ch_reference_segments)
        )


        // [sample_id, [final VCFs]] -> [sample_id, segment, vcf], the segment read back off
        // the {sample_id}__{ref}__{subtype}__{segment}_paf_final.vcf name post_process_vcf
        // writes. baseName, not simpleName: reference accessions carry dots.
        ch_final_vcf = post_process_vcf.out.vcf
            .transpose()
            .map { sample_id, vcf ->
                def segment = vcf.baseName.split('__')[-1].replaceFirst(/_paf_final$/, '')
                [sample_id, segment, vcf]
            }

        // Per-sample provenance. Inner joins: prepare_consensus, blastn and paftools_align
        // run errorStrategy 'ignore', so a sample failing any of them leaves the chain and
        // publishes no provenance file at all.
        ch_provenance = ch_provenance.combine(ch_pipeline_provenance).map{ it -> [it[0], [it[1]]] }
        ch_provenance = ch_provenance.join(hash_files.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
        ch_provenance = ch_provenance.join(prepare_consensus.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
        ch_provenance = ch_provenance.join(flumut.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
        ch_provenance = ch_provenance.join(blastn.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
        ch_provenance = ch_provenance.join(paftools_align.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
        ch_provenance = ch_provenance.join(paftools_call.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
        ch_provenance = ch_provenance.join(paftools_sort_vcf.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
        ch_provenance = ch_provenance.join(paftools_snpeff_annotation.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
        ch_provenance = ch_provenance.join(post_process_vcf.out.provenance).map{ it -> [it[0], it[1] << it[2]] }

        // --- Step 9: Scan the final VCF for a custom watchlist (optional) ---
        if (params.watchlist_dir != 'NO_FILE') {
            // checkIfExists so a mistyped --watchlist_dir fails loudly rather than
            // silently never scanning.
            ch_watchlists = Channel.fromPath(
                params.watchlist_dir + "/*watchlist.csv", checkIfExists: true
            )
            mutation_watch(ch_final_vcf.combine(ch_watchlists))

            // One task per (sample, segment, watchlist), so group back to one entry per
            // sample. remainder: true keeps a sample whose every scan failed, which an
            // inner join would drop along with the records it had already accumulated.
            ch_provenance = ch_provenance
                .join(mutation_watch.out.provenance.groupTuple(), remainder: true)
                .map{ it -> [it[0], (it[1] ?: []) + (it[2] ?: [])] }
        }

        collect_provenance(ch_provenance)

}
