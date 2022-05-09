// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// This step use GATK to count the number of varaint in each vcf file.
// The number of variant is one parameter used to filter samples based on quality.
// This step is followed by count_variants_gatk_2 to aggregate all the individual values into one file.

process count_variants_gatk {
        tag "${vcf}"

	input :
	file vcf
	file vcf_index 
	val assembly
	val batch
	val run
	
	output :
	file '*'

	script :
	"""
	gatk CountVariants \
	-V ${vcf} \
	-O ${vcf.simpleName}_gatk_count_variants
	"""
}

