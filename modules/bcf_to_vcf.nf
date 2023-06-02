// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// SNV Calling. 
// Split the multiallelic varaints (norm step) and transform the bcf into a vcf 
// Rename the varaints and compress the vcf into a vcf.gz
// Index the compressed vcf

process bcf_to_vcf {
	label 'conda_annotate'
        publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/SNV/", mode: 'copy'

	input :
	file bcf_file
	val assembly
	val batch
	val run

	output :
	path '*.vcf.gz', emit : vcf	

	script :
	"""
##	bcftools view ${bcf_file} | bgzip -c > ${bcf_file.simpleName}.vcf.gz
	bcftools norm -m -any -o ${bcf_file.simpleName}_norm.vcf ${bcf_file}
	bcftools annotate --set-id '%CHROM\\_%POS\\_%REF\\_%FIRST_ALT' -O z -o ${bcf_file.simpleName}.vcf.gz ${bcf_file.simpleName}_norm.vcf
	"""
}
