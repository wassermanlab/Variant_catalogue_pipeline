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
include { SV } from "./subworkflow/SV"

workflow{
	samples 	= Channel
		.fromFilePairs (params.reads)
		.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }    

	batch 		= params.batch
	assembly        = params.assembly
	run		= params.run
	outdir_ind 	= file (params.outdir_ind)

	main :
        ALN()
	QC_indiv(ALN.out.bam_sorted, ALN.out.bam_sorted_index)
        SNV(ALN.out.bam_sorted, ALN.out.bam_sorted_index, QC_indiv.out.mosdepth_output)
        MT(ALN.out.bam_sorted, ALN.out.bam_sorted_index)
	SV(ALN.out.bam_sorted, ALN.out.bam_sorted_index, SNV.out.sample_sex_file)
}

//Will have to add to the SV pipeline the output from SNV that defines the sex of each individual : SNV.out.{file_name}

// The main worflow can directly call the named workflows from the modules
//workflow {
//        ALN()
//	QC_indiv(ALN.out.bam_sorted, ALN.out.bam_sorted_index)
//        SNV(ALN.out.bam_sorted, ALN.out.bam_sorted_index, QC_indiv.out.mosdepth_output)
//        MT(ALN.out.bam_sorted, ALN.out.bam_sorted_index)
//}














