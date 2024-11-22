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
	file ref

	output :
	path '*_norm.vcf.gz', emit : vcf	
	path '*GLnexus_output.vcf.gz'

	script :
	"""
	# output an unmodified population vcf in compressed format
	bcftools view ${bcf_file} -Oz -o ${bcf_file.simpleName}_GLnexus_output.vcf.gz
	bcftools index -t ${bcf_file.simpleName}_GLnexus_output.vcf.gz

	# normalize/left align and split multi-allelic variants 
	bcftools norm -m -any -o ${bcf_file.simpleName}_norm_int.vcf -f ${ref} ${bcf_file}
	bcftools annotate --set-id '%CHROM\\_%POS\\_%REF\\_%FIRST_ALT' -O z -o ${bcf_file.simpleName}_norm.vcf.gz ${bcf_file.simpleName}_norm_int.vcf
	"""
}
