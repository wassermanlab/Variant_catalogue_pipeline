// Nextflow sub-workflow
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the sub-workflow goal and characteristics :
// Prepare the gnomad frequency files to extract only the necessary information

//This steps should be done only once and not re-run every time a sample is added to the BVL

// Load the modules for the ALN workflow

include { gnomad_frequency_table } from "./../modules/gnomad_frequency_table"

// Initialisation workflow
workflow Initialisation {

	// Load the parameters and files
	gnomad_SNV_vcf          = file (params.gnomad_SNV_vcf)
	gnomad_SNV_index        = file (params.gnomad_SNV_vcf_index)
	assembly        		= params.assembly
        chr                     = params.chrom


	main:
		gnomad_frequency_table(gnomad_SNV_vcf, gnomad_SNV_index, chr,assembly)
}


