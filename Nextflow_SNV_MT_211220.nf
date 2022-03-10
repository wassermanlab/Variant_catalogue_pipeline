#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// this prints the input parameters

log.info """
Dry-run for the Silent Genomes Project IBVL
=============================================
reads                           : ${params.reads}
reference                       : ${params.ref}
assembly                        : ${params.assembly}
run	                        : ${params.run}
batch	                        : ${params.batch}
"""

if (params.help) {
    log.info 'This is the dry-run for the IBVL Pipeline (Silent Genomes Project)'
    log.info 'Let\'s reduce health care disparities'
    exit 1
}


// Include the other workflow that themselves includes the modules
include { ALN } from "./subworkflow/ALN"
include { QC_indiv } from "./subworkflow/QC_indiv"
include { SNV } from "./subworkflow/SNV"
include { MT } from "./subworkflow/MT"


// The main worflow can directly call the named workflows from the modules
workflow {
        ALN()
	QC_indiv(ALN.out.bam_sorted, ALN.out.bam_sorted_index)
        SNV(ALN.out.bam_sorted, ALN.out.bam_sorted_index, QC_indiv.out.mosdepth_output)
        MT(ALN.out.bam_sorted, ALN.out.bam_sorted_index)
}














