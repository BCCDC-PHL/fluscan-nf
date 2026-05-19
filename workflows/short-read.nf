#!/usr/bin/env nextflow

/*
 *   Short-read Influenza A sequencing workflow
 *   -------------------------------------------
 *   Processes Illumina short-read FASTQ input through trimming, assembly,
 *   reference-based mapping, variant calling, and consensus generation.
 */

nextflow.enable.dsl = 2

include { short_read_hash_files }          from '../modules/short_read/hash_files.nf'
include { pipeline_provenance }            from '../modules/shared/provenance.nf'
include { collect_provenance }             from '../modules/shared/provenance.nf'
include { short_read_fastp }               from '../modules/short_read/qc.nf'
include { short_read_cutadapt }            from '../modules/short_read/qc.nf'
include { short_read_multiqc }             from '../modules/short_read/qc.nf'
include { short_read_fastqc }              from '../modules/short_read/qc.nf'
include { short_read_sample_reads }        from '../modules/short_read/assembly.nf'
include { short_read_assembly }            from '../modules/short_read/assembly.nf'
include { clade_calling }                  from '../modules/short_read/clade_calling.nf'
include { snp_calling }                    from '../modules/short_read/snp_calling.nf'
include { pull_genoflu }                   from '../modules/shared/genoflu.nf'
include { checkout_genoflu }               from '../modules/shared/genoflu.nf'
include { genoflu }                        from '../modules/shared/genoflu.nf'
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
include { short_read_summarize_vcf }       from '../modules/short_read/variant_calling.nf'
include { short_read_mutation_watch }      from '../modules/short_read/variant_calling.nf'
include { short_read_make_coverage_plot }  from '../modules/short_read/plotting.nf'


workflow SHORT_READ {

    def extract_segment = { sample_id, vcf ->
        def segment = vcf =~ /.+__([^_]+)__anno.vcf/
        println "DEBUG: extract_segment - sample_id: ${sample_id}, vcf: ${vcf}, segment: ${segment}"
        segment = segment[0] ? segment[0][1] : "unknown_segment"
        return [sample_id, segment, vcf]
    }

    ch_workflow_metadata = Channel.value([
        workflow.sessionId,
        workflow.runName,
        workflow.manifest.name,
        workflow.manifest.version,
        workflow.start,
    ])

    ch_pipeline_provenance = pipeline_provenance(ch_workflow_metadata)

    ch_primers = Channel.fromPath(params.primer_path)

    if (params.samplesheet_input != 'NO_FILE') {
        ch_fastq_input = Channel.fromPath(params.samplesheet_input).splitCsv(header: true).map{ it -> [it['ID'], it['R1'], it['R2']] }
    } else {
        ch_fastq_input = Channel.fromFilePairs( params.fastq_search_path, flat: true ).map{ it -> [it[0].split('_')[0], it[1], it[2]] }.unique{ it -> it[0] }
    }

    ch_aa_database = Channel.of([file(params.aa_database).parent, file(params.aa_database).name]).first()
    ch_nt_database = Channel.fromPath(params.nt_database)

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
        makeblastdb(ch_nt_database).first().set{ch_blastdb}
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

        if (params.snpeff_config != 'NO_FILE') {
            watchlist_files = Channel.fromPath(params.watchlist_dir + "/*watchlist.csv")
            short_read_snpeff_annotation(short_read_lofreq.out.vcf.combine(Channel.fromPath(params.snpeff_config)))

            ch_segment_vcfs = short_read_snpeff_annotation.out.vcfs.transpose().map(extract_segment)

            post_process_vcf(ch_segment_vcfs.join(blastn.out.ref))
            short_read_summarize_vcf(post_process_vcf.out.vcf)
            short_read_mutation_watch(short_read_summarize_vcf.out.tsv.combine(watchlist_files))
        }

        short_read_mask_low_coverage(short_read_sort_filter_sam.out.bam)
        short_read_make_consensus(short_read_lofreq.out.vcf.join(blastn.out.ref).join(short_read_mask_low_coverage.out.masked_bed))

        flumut(short_read_make_consensus.out.consensus)

        //Call clades for H1 and H3 samples
        //clade_calling(short_read_make_consensus.out.consensus)

        //genoflu(short_read_make_consensus.out.consensus)
        //snp_calling(short_read_make_consensus.out.consensus, ch_aa_database)

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
