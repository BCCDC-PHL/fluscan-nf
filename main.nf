#!/usr/bin/env nextflow

/*
 *   A nextflow wrapper for running FluScan
 *   -----------------------------------------
 */

nextflow.enable.dsl = 2

include { SHORT_READ } from './workflows/short-read.nf'
include { CONSENSUS }  from './workflows/consensus.nf'
include { NEXTCLADE }  from './workflows/nextclade.nf'


// prints to the screen and to the log
log.info """
  Flu Scanner Pipeline
  ===================================
  projectDir        : ${projectDir}
  launchDir         : ${launchDir}
  fastqInputDir     : ${params.fastq_input}
  fastaInputDir     : ${params.fasta_input}
  database          : ${params.nt_database}
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
        //CONSENSUS()
        NEXTCLADE()
    } else {
        error "ERROR: Please provide either --fastq_input (short-read workflow) or --fasta_input (consensus workflow)."
    }

}
