// Nextflow sub-workflow
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the sub-workflow goal and characteristics :
// Hail producing several graphs and filtering outliers samples and variants

// Load the modules for the SNV workflow

include { Hail_sample_QC } from "./../modules/Hail_sample_QC"
include { Hail_variant_QC } from "./../modules/Hail_variant_QC"


// Hail workflow

workflow Hail {

	// Load the parameters and files
	run             			= params.run
	batch           			= params.batch
	assembly        			= params.assembly
	assembly_hg				= params.assembly_hg

	// Workflow start here
	take : 
		SNV_vcf
	main :
		Hail_sample_QC(SNV_vcf, assembly, batch, run)
                Hail_variant_QC(Hail_sample_QC.out.vcf_sample_filtered, assembly, batch, run)
}

