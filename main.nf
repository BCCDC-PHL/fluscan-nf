#!/usr/bin/env nextflow

/*
 *   A nextflow wrapper for running FluScan
 *   -----------------------------------------
 */

include { SHORT_READ } from './workflows/short-read.nf'
include { FASTA_PREP } from './workflows/fasta-prep.nf'
include { PAFTOOLS }   from './workflows/paftools.nf'
include { NEXTCLADE }  from './workflows/nextclade.nf'
include { resolve_db } from './modules/shared/database.nf'
include { db_summary } from './modules/shared/database.nf'

// Resolved with no required components: this only checks that --db points at something
// that exists, so a typo fails before the banner. Each workflow re-resolves with the
// components it actually consumes.
def db = resolve_db(params.db, [])

// prints to the screen and to the log
log.info """
  Flu Scanner Pipeline
  ===================================
  projectDir        : ${projectDir}
  launchDir         : ${launchDir}
  fastqInputDir     : ${params.fastq_input}
  fastaInputDir     : ${params.fasta_input}
  samplesheet       : ${params.samplesheet_input}
  database          : ${db_summary(db)}
  primers           : ${params.primers}
  outdir            : ${params.outdir}
  pipeline run name : ${workflow.runName}
  pipeline version  : ${workflow.manifest.version}
  run_name          : ${params.run_name}
  user              : ${workflow.userName}
  Git repository    : ${workflow.repository}
  git commit id     : ${workflow.commitId}
  branch            : ${workflow.revision}
""".stripIndent()


workflow {

    if ([params.fastq_input, params.fasta_input, params.samplesheet_input].count { it != 'NO_FILE' } > 1) {
        error "ERROR: Provide only one of --fastq_input, --fasta_input or --samplesheet_input."
    }

    // A sample sheet's header row picks the workflow; samplesheet_channel validates the rest.
    def header = params.samplesheet_input != 'NO_FILE'
        ? (file(params.samplesheet_input, checkIfExists: true).splitCsv(strip: true)[0] ?: [])
        : []

    if (params.fastq_input != 'NO_FILE' || 'R1' in header) {
        SHORT_READ()
    } else if (params.fasta_input != 'NO_FILE' || 'FASTA' in header) {
        FASTA_PREP()
        PAFTOOLS(FASTA_PREP.out.cleaned, FASTA_PREP.out.provenance)
        NEXTCLADE(FASTA_PREP.out.cleaned, FASTA_PREP.out.provenance)
    } else if (params.samplesheet_input != 'NO_FILE') {
        error "ERROR: --samplesheet_input needs an ID,R1,R2 header (short-read workflow) or an ID,FASTA header (PAFTOOLS/NEXTCLADE workflow)."
    } else {
        error "ERROR: Please provide --fastq_input (short-read workflow), --fasta_input (PAFTOOLS/NEXTCLADE workflow) or --samplesheet_input (either)."
    }

}
