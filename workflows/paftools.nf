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
include { resolve_dbs }         from '../modules/shared/database.nf'
include { input_db_paths }      from '../modules/shared/sample_id.nf'

workflow PAFTOOLS {

    take:
        ch_cleaned              // [sample_id, cleaned consensus FASTA]
        ch_db                   // [sample_id, bundle name]
        ch_prep_provenance      // [sample_id, [hash, prepare, flumut]]

    main:

        // Reference selection is BLAST-driven and every called variant is SnpEff-annotated,
        // so both components are mandatory for this branch.
        def dbs = resolve_dbs(input_db_paths(), ['blastn', 'snpeff'])

        ch_workflow_metadata = Channel.fromList(dbs.values().collect { db -> [
            db.name,
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
        ] })

        ch_pipeline_provenance = pipeline_provenance(ch_workflow_metadata)    // [bundle name, record]

        makeblastdb(Channel.fromList(dbs.values().collect { db -> [db.name, db.blastn] }))

        // [sample_id, consensus, the sample's bundle's BLAST db]
        blastn(ch_cleaned.join(ch_db)
            .map { sample_id, consensus, db_name -> [db_name, sample_id, consensus] }
            .combine(makeblastdb.out.db, by: 0)
            .map { db_name, sample_id, consensus, blastdb -> [sample_id, consensus, blastdb] })

        // [sample_id, consensus, BLAST references, snpeff/]. The whole snpeff/ directory is
        // staged so the config and the database it names travel together; the process
        // overrides data.dir with -dataDir.
        paftools_variants(ch_cleaned.join(blastn.out.ref).join(ch_db)
            .map { sample_id, consensus, reference, db_name -> [sample_id, consensus, reference, dbs[db_name].snpeff_dir] })

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
        ch_provenance = ch_prep_provenance.join(ch_db)
            .map { sample_id, records, db_name -> [db_name, sample_id, records] }
            .combine(ch_pipeline_provenance, by: 0)
            .map { db_name, sample_id, records, pipeline -> [sample_id, records, pipeline] }
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
