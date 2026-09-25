#!/usr/bin/env nextflow

/*
 *   Short-read Influenza A sequencing workflow
 *   -------------------------------------------
 *   Processes Illumina short-read FASTQ input through trimming, assembly,
 *   reference-based mapping, variant calling, and consensus generation.
 */

include { hash_files as short_read_hash_files } from '../modules/shared/hash_files.nf'
include { pipeline_provenance }            from '../modules/shared/provenance.nf'
include { collect_provenance }             from '../modules/shared/provenance.nf'
include { short_read_fastp }               from '../modules/short_read/qc.nf'
include { short_read_cutadapt }            from '../modules/short_read/qc.nf'
include { short_read_multiqc }             from '../modules/short_read/qc.nf'
include { short_read_fastqc }              from '../modules/short_read/qc.nf'
include { short_read_sample_reads }        from '../modules/short_read/assembly.nf'
include { short_read_assembly }            from '../modules/short_read/assembly.nf'
include { flumut }                         from '../modules/shared/flumut.nf'

include { makeblastdb }                    from '../modules/shared/blast.nf'
include { blastn }                         from '../modules/shared/blast.nf'
include { short_read_map_reads_bwa }       from '../modules/short_read/alignment.nf'
include { short_read_map_reads_minimap2 }  from '../modules/short_read/alignment.nf'
include { short_read_sort_filter_sam }     from '../modules/short_read/alignment.nf'
include { short_read_make_pileup }         from '../modules/short_read/alignment.nf'
include { short_read_mask_low_coverage }   from '../modules/short_read/alignment.nf'
include { short_read_make_consensus }      from '../modules/short_read/alignment.nf'
include { short_read_lofreq }              from '../modules/short_read/variant_calling.nf'
include { short_read_snpeff_annotation }   from '../modules/short_read/variant_calling.nf'
include { post_process_vcf }               from '../modules/shared/post_process.nf'
include { split_fasta as split_fasta_reference } from '../modules/shared/consensus_qc.nf'
include { mutation_watch }                 from '../modules/shared/mutation_watch.nf'
include { short_read_make_coverage_plot }  from '../modules/short_read/plotting.nf'
include { resolve_db }                     from '../modules/shared/database.nf'
include { samplesheet_channel }            from '../modules/shared/sample_id.nf'


workflow SHORT_READ {

    // Only the BLAST reference is required here. SnpEff annotation (and the MNV
    // correction and watchlist scan that depend on its output) runs when the bundle
    // carries a snpeff/ component.
    def db = resolve_db(params.db, ['blastn'])

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

    ch_primers = Channel.fromPath(params.primers, checkIfExists: true)

    if (params.samplesheet_input != 'NO_FILE') {
        ch_fastq_input = samplesheet_channel(params.samplesheet_input, ['R1', 'R2'])
    } else {
        ch_fastq_input = Channel.fromFilePairs( params.fastq_search_path, flat: true ).map{ it -> [it[0].split('_')[0], it[1], it[2]] }.unique{ it -> it[0] }
    }

    ch_blastn_db = Channel.fromPath(db.blastn)

    main:

        // Provenance channel starts with just the sample IDs
        // These will be joined to various provenance files as they are generated
        ch_provenance = ch_fastq_input.map{ it -> it[0] }

        // Generate hashes for input files
        short_read_hash_files(ch_fastq_input.map{ it -> [it[0], [it[1], it[2]]] }.combine(Channel.of("fastq_input")))

        // Clean up reads - remove adapters (fastp) and primers (cutadapt)
        short_read_fastp(ch_fastq_input)
        short_read_cutadapt(short_read_fastp.out.trimmed_reads.combine(ch_primers))
        short_read_fastqc(short_read_cutadapt.out.primer_trimmed_reads)

        short_read_sample_reads(short_read_cutadapt.out.primer_trimmed_reads)

        // Run assembly
        short_read_assembly(short_read_sample_reads.out.reads)

        // Collect all the relevant files for multiqc
        ch_fastqc_collected = short_read_fastqc.out.zip.map{ it -> [it[1], it[2]]}.collect()
        short_read_multiqc(short_read_fastp.out.json.mix( short_read_cutadapt.out.log, ch_fastqc_collected ).collect().ifEmpty([]) )

        // Build BLAST database, then search assembled contigs against it to select references
        makeblastdb(ch_blastn_db).first().set{ch_blastdb}
        blastn(short_read_assembly.out.contigs, ch_blastdb)

        // NOTE: Currently mapping raw reads (ch_fastq_input). Consider using
        // short_read_cutadapt.out.primer_trimmed_reads instead to avoid adapter/primer artifacts.
        if (params.aligner == 'minimap2') {
            short_read_map_reads_minimap2(ch_fastq_input.join(blastn.out.ref))
            ch_sam = short_read_map_reads_minimap2.out.sam
        } else {
            short_read_map_reads_bwa(ch_fastq_input.join(blastn.out.ref))
            ch_sam = short_read_map_reads_bwa.out.sam
        }
        short_read_sort_filter_sam(ch_sam)
        short_read_make_pileup(short_read_sort_filter_sam.out.bam)
        short_read_make_coverage_plot(short_read_sort_filter_sam.out.depth)
        short_read_lofreq(short_read_sort_filter_sam.out.bam.join(blastn.out.ref))

        if (db.has_snpeff) {
            // The whole snpeff/ directory is staged so the config and the database it
            // names travel together; the process overrides data.dir with -dataDir.
            short_read_snpeff_annotation(short_read_lofreq.out.vcf.combine(Channel.fromPath(db.snpeff_dir)))

            // Split each sample's multi-segment BLAST reference into per-segment FASTAs
            // named {sample_id}_ref_{segment}.s.fasta — the lookup post_process_vcf expects
            ch_reference_segments = split_fasta_reference(blastn.out.ref)


            // Fully aggregated per sample: [sample_id, [segment VCFs], [per-segment refs]]
            ch_segment_vcfs = short_read_snpeff_annotation.out.vcfs.join(ch_reference_segments)

            post_process_vcf(ch_segment_vcfs)

            // Watchlist scanning reads the final VCF, so it keys only on a watchlist being
            // set -- but the final VCF only exists inside this snpeff block, which is what
            // keeps the call here.
            if (params.watchlist_dir != 'NO_FILE') {
                // checkIfExists so a mistyped --watchlist_dir fails loudly. Built inside the
                // guard: an unguarded fromPath('NO_FILE/*watchlist.csv') yields an empty
                // channel, and the scan then silently does not run at all.
                ch_watchlists = Channel.fromPath(
                    params.watchlist_dir + "/*watchlist.csv", checkIfExists: true
                )
                // [sample_id, vcf] -> [sample_id, segment, vcf]; the segment is carried as
                // '' because short-read post_process_vcf emits one VCF per sample, not per
                // segment, and an invented label would be reported in the output.
                ch_watch_in = post_process_vcf.out.vcf
                    .transpose()
                    .map { sample_id, vcf -> [sample_id, '', vcf] }
                    .combine(ch_watchlists)
                mutation_watch(ch_watch_in)
            }
        }

        short_read_mask_low_coverage(short_read_sort_filter_sam.out.bam)
        short_read_make_consensus(short_read_lofreq.out.vcf.join(blastn.out.ref).join(short_read_mask_low_coverage.out.masked_bed))

        flumut(short_read_make_consensus.out.consensus)

        // Provenance collection processes
        // The basic idea is to build up a channel with the following structure:
        // [sample_id, [provenance_file_1.yml, provenance_file_2.yml, provenance_file_3.yml...]]
        // ...and then concatenate them all together in the 'collect_provenance' process.
        ch_provenance = ch_provenance.combine(ch_pipeline_provenance).map{ it ->      [it[0], [it[1]]] }
        ch_provenance = ch_provenance.join(short_read_hash_files.out.provenance).map{ it ->      [it[0], it[1] << it[2]] }
        ch_provenance = ch_provenance.join(short_read_fastp.out.provenance).map{ it ->           [it[0], it[1] << it[2]] }
        ch_provenance = ch_provenance.join(short_read_cutadapt.out.provenance).map{ it ->        [it[0], it[1] << it[2]] }
        ch_provenance = ch_provenance.join(short_read_assembly.out.provenance).map{ it ->        [it[0], it[1] << it[2]] }
        ch_provenance = ch_provenance.join(flumut.out.provenance).map{ it ->         [it[0], it[1] << it[2]] }
        collect_provenance(ch_provenance)

}
