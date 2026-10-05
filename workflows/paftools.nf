#!/usr/bin/env nextflow

/*
 *   PAFTOOLS workflow (consensus FASTA input)
 *   -----------------------------------------
 *   Takes the consensus FASTAs FASTA_PREP has cleaned through subtype-aware reference
 *   selection via BLAST, minimap2/paftools variant calling, bcftools VCF sorting, SnpEff
 *   annotation, and summary extraction.
 *
 *   Reference database header format (<--db>/blastn_db.fa):
 *     >accession|subtype|segment
 *   e.g.
 *     >NC_026438.1|H1N1|PB2
 *
 */

include { pipeline_provenance }                          from '../modules/shared/provenance.nf'
include { collect_provenance }                           from '../modules/shared/provenance.nf'
include { split_fasta as split_fasta_reference }         from '../modules/shared/consensus_qc.nf'
include { split_fasta as split_fasta_consensus }         from '../modules/shared/consensus_qc.nf'
include { paftools_align }                               from '../modules/paftools/alignment.nf'
include { consensus_coverage_qc }                        from '../modules/paftools/alignment.nf'
include { makeblastdb }                                  from '../modules/shared/blast.nf'
include { blastn }                                       from '../modules/shared/blast.nf'
include { paftools_call }                                from '../modules/paftools/variant_calling.nf'
include { paftools_sort_vcf }                            from '../modules/paftools/variant_calling.nf'
include { paftools_snpeff_annotation }                   from '../modules/paftools/variant_calling.nf'
include { paftools_to_freebayes_vcf }                    from '../modules/paftools/variant_calling.nf'
include { post_process_vcf }                             from '../modules/shared/post_process.nf'
include { mutation_watch }                               from '../modules/shared/mutation_watch.nf'
include { resolve_db }                                   from '../modules/shared/database.nf'



workflow PAFTOOLS {

    take:
        ch_cleaned              // [sample_id, cleaned consensus FASTA]
        ch_prep_provenance      // [sample_id, [hash, prepare, flumut]]

    main:

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

        ch_blastn_db = Channel.fromPath(db.blastn)

        // The whole snpeff/ directory is staged so the config and the database it names
        // travel together; the process overrides data.dir with -dataDir.
        ch_snpeff_dir = Channel.fromPath(db.snpeff_dir, checkIfExists: true)

        // --- Step 2: Build BLAST database (once, shared across all samples) ---
        makeblastdb(ch_blastn_db).first().set{ch_blastdb}

        // --- Step 3: BLAST cleaned consensus against reference database ---
        blastn(ch_cleaned, ch_blastdb)

        ch_consensus_segments = split_fasta_consensus(ch_cleaned)
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


        // [sample_id, [final VCFs]] -> [sample_id, [segments], [vcfs]], each segment read back
        // off the {sample_id}__{ref}__{subtype}__{segment}_paf_final.vcf name post_process_vcf
        // writes. baseName, not simpleName: reference accessions carry dots.
        ch_final_vcfs = post_process_vcf.out.vcf
            .map { sample_id, vcfs ->
                vcfs = [vcfs].flatten()
                def segments = vcfs.collect { it.baseName.split('__')[-1].replaceFirst(/_paf_final$/, '') }
                [sample_id, segments, vcfs]
            }

        // --- Step 9: Scan the final VCFs for custom watchlists (optional) ---
        if (params.watchlist_dir != 'NO_FILE') {
            // checkIfExists so a mistyped --watchlist_dir fails loudly rather than
            // silently never scanning.
            ch_watchlists = Channel.fromPath(
                params.watchlist_dir + "/*watchlist.csv", checkIfExists: true
            )
            mutation_watch(ch_final_vcfs, ch_watchlists.collect())
        }

        // Per-sample provenance. Inner joins: prepare_consensus, blastn and paftools_align
        // run errorStrategy 'ignore', so a sample failing any of them leaves the chain and
        // publishes no provenance file at all.
        ch_provenance = ch_prep_provenance
            .combine(ch_pipeline_provenance)
            .map{ sample_id, records, pipeline -> [sample_id, [pipeline] + records] }
        ch_provenance = ch_provenance.join(blastn.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
        ch_provenance = ch_provenance.join(paftools_align.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
        ch_provenance = ch_provenance.join(paftools_call.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
        ch_provenance = ch_provenance.join(paftools_sort_vcf.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
        ch_provenance = ch_provenance.join(paftools_snpeff_annotation.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
        ch_provenance = ch_provenance.join(post_process_vcf.out.provenance).map{ it -> [it[0], it[1] << it[2]] }

        if (params.watchlist_dir != 'NO_FILE') {
            // One scan task per sample. remainder: true keeps samples whose scan failed.
            ch_provenance = ch_provenance
                .join(mutation_watch.out.provenance, remainder: true)
                .map{ sample_id, records, watch -> [sample_id, (records ?: []) + (watch ? [watch] : [])] }
        }

        collect_provenance(ch_provenance)

}
