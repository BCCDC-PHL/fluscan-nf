#!/usr/bin/env nextflow

/*
 *   A nextflow wrapper for running FluScan
 *   -----------------------------------------
 */

nextflow.enable.dsl = 2

include { SHORT_READ } from './workflows/short-read.nf'
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

    if (params.fastq_input != 'NO_FILE') {
        SHORT_READ()
    } else if (params.fasta_input != 'NO_FILE') {
        PAFTOOLS()
        NEXTCLADE()
    } else {
        error "ERROR: Please provide either --fastq_input (short-read workflow) or --fasta_input (PAFTOOLS/NEXTCLADE workflow)."
    }

}
