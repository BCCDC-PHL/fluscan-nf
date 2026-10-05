#!/usr/bin/env nextflow

/*
 *   PAFTOOLS workflow (consensus FASTA input)
 *   -----------------------------------------
 *   Takes the consensus FASTAs FASTA_PREP has cleaned through subtype-aware reference
 *   selection via BLAST, then minimap2/paftools variant calling, SnpEff annotation and
 *   MNV codon correction in one task per sample.
 *
 *   Reference database header format (<--db>/blastn_db.fa):
 *     >accession|subtype|segment
 *   e.g.
 *     >NC_026438.1|H1N1|PB2
 *
 */

include { pipeline_provenance } from '../modules/shared/provenance.nf'
include { collect_provenance }  from '../modules/shared/provenance.nf'
include { makeblastdb }         from '../modules/shared/blast.nf'
include { blastn }              from '../modules/shared/blast.nf'
include { paftools_variants }   from '../modules/paftools/paftools.nf'
include { mutation_watch }      from '../modules/shared/mutation_watch.nf'
include { resolve_db }          from '../modules/shared/database.nf'

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

        // The whole snpeff/ directory is staged so the config and the database it names
        // travel together; the process overrides data.dir with -dataDir.
        ch_snpeff_dir = Channel.fromPath(db.snpeff_dir, checkIfExists: true)

        makeblastdb(Channel.fromPath(db.blastn)).first().set{ch_blastdb}

        blastn(ch_cleaned, ch_blastdb)

        // [sample_id, consensus, BLAST references, snpeff/]
        paftools_variants(ch_cleaned.join(blastn.out.ref).combine(ch_snpeff_dir))

        paftools_variants.out.coverage
            .map{ it -> it[1] }
            .flatten()                              // one Path per emission -> avoids the [name, content...] branch
            .collectFile(
                name: 'paftools_qc_summary.tsv',
                storeDir: "${params.outdir}/aggregate",
                keepHeader: true,                   // keep header from first file only
                skip: 1,                            // drop the header line from the rest
                sort: false
            )

        // [sample_id, [final VCFs]] -> [sample_id, [segments], [vcfs]], each segment read back
        // off the {sample_id}__{ref}__{subtype}__{segment}_paf_final.vcf name.
        // baseName, not simpleName: reference accessions carry dots.
        ch_final_vcfs = paftools_variants.out.vcf
            .map { sample_id, vcfs ->
                vcfs = [vcfs].flatten()
                def segments = vcfs.collect { it.baseName.split('__')[-1].replaceFirst(/_paf_final$/, '') }
                [sample_id, segments, vcfs]
            }

        if (params.watchlist_dir != 'NO_FILE') {
            // checkIfExists so a mistyped --watchlist_dir fails loudly rather than
            // silently never scanning.
            ch_watchlists = Channel.fromPath(
                params.watchlist_dir + "/*watchlist.csv", checkIfExists: true
            )
            mutation_watch(ch_final_vcfs, ch_watchlists.collect())
        }

        // Per-sample provenance. Inner joins: prepare_consensus and blastn run errorStrategy
        // 'ignore', so a sample failing either leaves the chain and publishes no provenance.
        ch_provenance = ch_prep_provenance
            .combine(ch_pipeline_provenance)
            .join(blastn.out.provenance)
            .join(paftools_variants.out.provenance)
            .map{ sample_id, records, pipeline, blast, paftools -> [sample_id, [pipeline] + records + [blast, paftools]] }

        if (params.watchlist_dir != 'NO_FILE') {
            // One scan task per sample. remainder: true keeps samples whose scan failed.
            ch_provenance = ch_provenance
                .join(mutation_watch.out.provenance, remainder: true)
                .map{ sample_id, records, watch -> [sample_id, (records ?: []) + (watch ? [watch] : [])] }
        }

        collect_provenance(ch_provenance)

}
