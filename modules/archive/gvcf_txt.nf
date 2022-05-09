// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// SNV Calling. Create a list with all the individual vcf that are included in the next step of joint calling 
// (!) Need to be done once all calls have been performed

process gvcfs_txt {

        publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/GLnexus/", mode: 'copy'

	input :
	path deepvariant_gvcf
	val assembly
	val batch
	val run

	output :
	path 'gvcfs.txt'

	script:
	"""
	find $params.outdir_ind/${assembly}/${batch}/${run}/DeepVariant/ -name "*.g.vcf.gz" > gvcfs.txt
	"""
}
