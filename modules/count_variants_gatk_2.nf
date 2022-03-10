// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// This step follow count_variants_gatk
// count_variants_gatk gives one output by sample, this steps aggregates all the outputs into one file, that is then laoded into R to create graphs for sample QC

process count_variants_gatk_2 {

	publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/QC/Out_MultiQC/variant_count/", mode: 'copy'

	input :
	path gatk_count_variant
	val assembly
	val batch
	val run
	
	output :
	file 'SNP_count_GATK.tsv'

	script :
	"""
	#For the number of snps from eahc vcf outputted by gatk, need to put them all in one file
        grep '.*' *_gatk_count_variants > SNP_count_GATK.tsv
	"""
}

