#!/usr/bin/env nextflow

/*
 *   Short-read Influenza A sequencing workflow
 *   -------------------------------------------
 *   Processes Illumina short-read FASTQ input through trimming, assembly,
 *   reference-based mapping, variant calling, and consensus generation.
 */

nextflow.enable.dsl = 2

include { hash_files }          from '../modules/hash_files.nf'
include { pipeline_provenance } from '../modules/provenance.nf'
include { collect_provenance }  from '../modules/provenance.nf'
include { fastp }               from '../modules/short_read_qc.nf'
include { cutadapt }            from '../modules/short_read_qc.nf'
include { multiqc }             from '../modules/short_read_qc.nf'
include { fastqc }              from '../modules/short_read_qc.nf'
include { sample_reads }        from '../modules/assembly.nf'
include { assembly }            from '../modules/assembly.nf'
include { clade_calling }       from '../modules/clade_calling.nf'
include { snp_calling }         from '../modules/snp_calling.nf'
include { pull_genoflu }        from '../modules/genoflu.nf'
include { checkout_genoflu }    from '../modules/genoflu.nf'
include { genoflu }             from '../modules/genoflu.nf'

include { makeblastdb }         from '../modules/blast.nf'
include { blastn }              from '../modules/blast.nf'
include { map_reads }           from '../modules/alignment.nf'
include { sort_filter_sam }     from '../modules/alignment.nf'
include { make_pileup }         from '../modules/alignment.nf'
include { mask_low_coverage }   from '../modules/alignment.nf'
include { make_consensus }      from '../modules/alignment.nf'
include { freebayes }           from '../modules/variant_calling.nf'
include { snpeff_annotation }   from '../modules/variant_calling.nf'
include { reformat_vcf }        from '../modules/variant_calling.nf'
include { mutation_watch }      from '../modules/variant_calling.nf'
include { make_coverage_plot }  from '../modules/plotting.nf'


workflow SHORT_READ {

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
        hash_files(ch_fastq_input.map{ it -> [it[0], [it[1], it[2]]] }.combine(Channel.of("fastq_input")))

        // Clean up reads - remove adapters (fastp) and primers (cutadapt)
        fastp(ch_fastq_input)
        cutadapt(fastp.out.trimmed_reads.combine(ch_primers))
        fastqc(cutadapt.out.primer_trimmed_reads)

        sample_reads(cutadapt.out.primer_trimmed_reads)

        // Run assembly 
        assembly(sample_reads.out.reads)

        // Collect all the relevant files for multiqc
        ch_fastqc_collected = fastqc.out.zip.map{ it -> [it[1], it[2]]}.collect()
        multiqc(fastp.out.json.mix( cutadapt.out.log, ch_fastqc_collected ).collect().ifEmpty([]) )

        // Build BLAST database, then search assembled contigs against it to select references
        makeblastdb(ch_nt_database).first().set{ch_blastdb}
        blastn(assembly.out.contigs, ch_blastdb)

        // NOTE: Currently mapping raw reads (ch_fastq_input). Consider using
        // cutadapt.out.primer_trimmed_reads instead to avoid adapter/primer artifacts.
        map_reads(ch_fastq_input.join(blastn.out.ref))
        sort_filter_sam(map_reads.out.sam)
        make_pileup(sort_filter_sam.out.bam)
        make_coverage_plot(sort_filter_sam.out.depth)
        freebayes(sort_filter_sam.out.bam.join(blastn.out.ref))

        if (params.snpeff_config != 'NO_FILE') {
            watchlist_files = Channel.fromPath(params.watchlist_dir + "/*watchlist.csv")
            snpeff_annotation(freebayes.out.vcf.combine(Channel.fromPath(params.snpeff_config)))
            reformat_vcf(snpeff_annotation.out.vcf)
            mutation_watch(reformat_vcf.out.tsv.combine(watchlist_files))
        }

        mask_low_coverage(sort_filter_sam.out.bam)
        make_consensus(freebayes.out.vcf.join(blastn.out.ref).join(mask_low_coverage.out.masked_bed))

        //Call clades for H1 and H3 samples
        //clade_calling(make_consensus.out.consensus)

        //genoflu(make_consensus.out.consensus)
        //snp_calling(make_consensus.out.consensus, ch_aa_database)

        // Provenance collection processes
        // The basic idea is to build up a channel with the following structure:
        // [sample_id, [provenance_file_1.yml, provenance_file_2.yml, provenance_file_3.yml...]]
        // ...and then concatenate them all together in the 'collect_provenance' process.
        ch_provenance = ch_provenance.combine(ch_pipeline_provenance).map{ it ->      [it[0], [it[1]]] }
        ch_provenance = ch_provenance.join(hash_files.out.provenance).map{ it ->      [it[0], it[1] << it[2]] }
        ch_provenance = ch_provenance.join(fastp.out.provenance).map{ it ->           [it[0], it[1] << it[2]] }
        ch_provenance = ch_provenance.join(cutadapt.out.provenance).map{ it ->        [it[0], it[1] << it[2]] }
        ch_provenance = ch_provenance.join(assembly.out.provenance).map{ it ->        [it[0], it[1] << it[2]] }
        collect_provenance(ch_provenance)

}
