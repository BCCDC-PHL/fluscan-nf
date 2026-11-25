#!/usr/bin/env nextflow

/*
 *   Consensus FASTA preparation, shared by PAFTOOLS and NEXTCLADE
 *   -------------------------------------------------------------
 *   Hashes each input FASTA, validates and cleans it, aggregates the consensus QC rows and
 *   runs FluMut, once per sample however many FASTA branches consume the result. With
 *   --db_pool, a sample without a sample sheet DB cell gets the pool bundle its HA and NA
 *   match best (select_db); a sample with no confident match goes to neither branch.
 *
 *   The sample name comes from the input filename (modules/shared/sample_id.nf); input
 *   headers supply only the segment label. Accepted header formats:
 *     Underscore-delimited (legacy): >SAMPLE_ID_segmentN_SEGMENT
 *       e.g. >WIN-AH-2024-FAV-0415-1-OS_segment1_PB2
 *     Pipe-delimited (GISAID-style): >EPI_ISL_ACCESSION|subtype|SEGMENT|...|clade_SEGMENT
 *       e.g. >EPI_ISL_20399898|A_/_H5N1|NA|6|EPI5261179|2.3.4.4b_NA
 *   Cleaned headers are normalised to {sample_id}_{SEGMENT} by prepare_consensus.
 */

include { hash_files }              from '../modules/shared/hash_files.nf'
include { prepare_consensus }       from '../modules/shared/consensus_qc.nf'
include { flumut }                  from '../modules/shared/flumut.nf'
include { fasta_sample_channel }    from '../modules/shared/sample_id.nf'
include { input_db_paths }          from '../modules/shared/sample_id.nf'
include { resolve_dbs }             from '../modules/shared/database.nf'
include { resolve_pool }            from '../modules/shared/database.nf'
include { makeblastdb }             from '../modules/shared/blast.nf'
include { select_db }               from '../modules/shared/select_db.nf'

workflow FASTA_PREP {

    main:

        // Consensus FASTAs from the sample sheet or the input directory, with each sample's
        // bundle. Each file is expected to contain 8 segments for a single sample.
        def dbs = resolve_dbs(input_db_paths(), [])
        ch_samples = fasta_sample_channel()      // [sample_id, FASTA, bundle path or null, input ID]
        ch_fasta_input = ch_samples.map { sample_id, fasta, db_path, input_id -> [sample_id, fasta] }
        ch_named_db = ch_samples.filter { it[2] }.map { sample_id, fasta, db_path, input_id ->
            [sample_id, dbs.values().find { it.root == file(db_path) }.name]
        }

        // Hash the input FASTA. One file per sample, so no list-wrapping here.
        hash_files(ch_fasta_input.combine(Channel.of("fasta_input")))

        // Validate and clean input consensus FASTA.
        prepare_consensus(ch_fasta_input)

        // One run-level table of the per-segment consensus QC rows.
        prepare_consensus.out.qc
            .map{ it -> it[1] }
            .collectFile(
                name: 'consensus_qc_summary.tsv',
                storeDir: "${params.outdir}/aggregate",
                keepHeader: true,                   // keep header from first file only
                skip: 1,                            // drop the header line from the rest
                sort: false
            )

        // HPAI mutation screening (per sample).
        flumut(prepare_consensus.out.cleaned)

        // Inner joins: prepare_consensus runs errorStrategy 'ignore', so a sample it rejects
        // leaves every branch's provenance chain.
        ch_provenance = hash_files.out.provenance
            .join(prepare_consensus.out.provenance)
            .join(flumut.out.provenance)
            .map { sample_id, hash, prepare, mutations -> [sample_id, [hash, prepare, mutations]] }

        ch_db = ch_named_db
        if (params.db_pool != 'NO_FILE') {
            // One BLAST db over every pool bundle, headers prefixed bundle|accession|subtype|segment.
            def pool = resolve_pool(params.db_pool, [])
            ch_pool_fasta = Channel.fromList(pool.values() as List)
                .collectFile(name: 'db_pool.fa', sort: true) { db -> db.blastn.text.trim().replaceAll(/(?m)^>/, ">${db.name}|") + '\n' }
            makeblastdb(ch_pool_fasta.map { fasta -> ['db_pool', fasta] })

            // Cleaned samples without a named bundle -> [sample_id, consensus, pool BLAST db]
            select_db(prepare_consensus.out.cleaned
                .join(ch_samples.filter { !it[2] }.map { [it[0], 'pool'] })
                .combine(makeblastdb.out.db)
                .map { sample_id, consensus, pool_label, db_name, blastdb -> [sample_id, consensus, blastdb] })

            select_db.out.selection
                .map { it[1] }
                .collectFile(
                    name: 'db_selection.tsv',
                    storeDir: "${params.outdir}/aggregate",
                    keepHeader: true,
                    skip: 1,
                    sort: true
                )

            // An empty bundle file is an unassigned sample: it leaves both branches here.
            ch_db = ch_named_db.mix(select_db.out.db
                .map { sample_id, db_file -> [sample_id, db_file.text.trim()] }
                .filter { it[1] })

            ch_provenance = ch_provenance
                .join(select_db.out.provenance, remainder: true)
                .filter { it[1] != null }
                .map { sample_id, records, selection -> [sample_id, records + (selection ? [selection] : [])] }
        }

        // Input ID to sample_id and bundle, for joining results back to the submitter's IDs.
        ch_samples.join(ch_db)
            .map { sample_id, fasta, db_path, input_id, db_name -> "${input_id},${sample_id},${db_name}\n" }
            .collectFile(
                name: 'sample_ids.csv',
                storeDir: "${params.outdir}/aggregate",
                seed: "input_id,sample_id,db\n",
                sort: true
            )

    emit:
        cleaned    = prepare_consensus.out.cleaned    // [sample_id, cleaned consensus FASTA]
        db         = ch_db                            // [sample_id, bundle name]
        provenance = ch_provenance                    // [sample_id, [hash, prepare, flumut(, select_db)]]
}
