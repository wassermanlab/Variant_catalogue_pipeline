#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// this prints the input parameters

log.info """
Varaint catalogue pipeline - Nextflow
=============================================
reads                           : ${params.reads}
reference                       : ${params.ref}
assembly                        : ${params.assembly}
run	                        : ${params.run}
batch	                        : ${params.batch}
"""

if (params.help) {
    log.info """
    Usage:
        This is the variant catalogue pipeline
	This pipeline was developped to generate the IBVL (Silent Genomes Project)
	Let\'s reduce health care disparities
	
        The typical command for running the pipeline is as follows:
        nextflow run wassermanlab/Variant_catalogue_pipeline -r main -profile GRCh38

        Mandatory arguments:
         -profile                      Allows to choose which reference genome to use, in this version, only GRCh37 and GRCh38 are available.

       Optional arguments:
        --help                         This usage statement
    """
    exit 1
}


// Include the other workflow that themselves includes the modules
include { Initialisation } from "./subworkflow/Initialisation"
include { Mapping } from  "./subworkflow/Mapping"
include { SNV } from "./subworkflow/SNV"
include { SV } from "./subworkflow/SV"
include { MT } from "./subworkflow/MT"

workflow{
	samples 	= Channel
		.fromFilePairs (params.reads)
		.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }    

	batch 		= params.batch
	assembly        = params.assembly
	run		= params.run
	outdir_ind 	= file (params.outdir_ind)

	main :
	//Initialisation()
        Mapping()
        SNV(Mapping.out.bam_sorted, Mapping.out.bam_sorted_index)
        SV(Mapping.out.bam_sorted, Mapping.out.bam_sorted_index, SNV.out.sample_sex_file)
	MT(Mapping.out.bam_sorted, Mapping.out.bam_sorted_index, Mapping.out.mosdepth_output)
}
