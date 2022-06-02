// Nextflow sub-workflow
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the sub-workflow goal and characteristics :
// Call the SNV variants
// Include some quality controls (QC) steps
//	- Plink which defines the sex of each sample based on seevral variables

// Load the modules for the SNV workflow

include { deepvariant_call } from "./../modules/deepvariant"
include { list_vcfs_txt } from "./../modules/list_vcfs_txt"
include { GLnexus_cli } from "./../modules/GLnexus"
include { bcf_to_vcf } from "./../modules/bcf_to_vcf"

include { plink_sex_inference } from "./../modules/plink_sex_inference"
include { sample_QC } from "./../modules/sample_QC"

// SNV workflow

workflow SNV {

	// Load the parameters and files
	run             			= params.run
	batch           			= params.batch
	assembly        			= params.assembly
	assembly_hg				= params.assembly_hg
	reference       			= file (params.ref)
	reference_index 			= file (params.ref_index)
        SNV                                     = params.SNV

	// Workflow start here
	take : 
		bam
		bai
		mosdepth

	main :
		// Sample specific (Do not need to be run for a previously processed sample)
		deepvariant_call(reference, reference_index, bam, bai, assembly, batch, run)

		// Aggregated steps (Need to be run everytime a new sample is added to the cohort)
		list_vcfs_txt(deepvariant_call.out.deepvariant_gvcf.collect(), assembly, batch, run, SNV)
		GLnexus_cli(list_vcfs_txt.out, run)
		bcf_to_vcf(GLnexus_cli.out, assembly, batch, run)

                plink_sex_inference(bcf_to_vcf.out.vcf, assembly_hg, assembly, batch, run)
		sample_QC(plink_sex_inference.out, assembly, batch, run, mosdepth)

	emit :
		sample_sex_file=sample_QC.out.sample_QC_file
		SNV_vcf = bcf_to_vcf.out.vcf
}
