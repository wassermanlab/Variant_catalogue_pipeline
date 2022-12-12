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
    log.info 'Pipeline designed to generate a genetic variant catalogue'
    log.info 'Pipeline developped to generate the IBVL (Silent Genomes Project)'
    log.info 'Let\'s reduce health care disparities'
    exit 1
}


// Include the other workflow that themselves includes the modules
//include { Initialisation } from "./subworkflow/Initialisation"
include { Mapping } from  "./subworkflow/Mapping"
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
//	Initialisation()
        Mapping()
	SNV(Mapping.out.bam_sorted, Mapping.out.bam_sorted_index)
        MT(Mapping.out.bam_sorted, Mapping.out.bam_sorted_index, Mapping.out.mosdepth_output)
	SV(Mapping.out.bam_sorted, Mapping.out.bam_sorted_index, SNV.out.sample_sex_file)
}
