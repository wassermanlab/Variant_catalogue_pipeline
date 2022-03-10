// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// Use plink to define the F coefficient estimates (i.e. (<observed hom. count> - <expected count>) / (<total observations> - <expected count>))
// The first step split off the X chromosome pseudo-autosomal region
// The second step calculate the F coefficient for each sample
// This value, as well as the X and Y normalized coverage are used to define the sample sex

process plink_sex_inference {
        tag "$vcf"

	publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/QC/Out_MultiQC/", mode: 'copy'

	input :
	path vcf
	val assembly
	val batch
	val run
	
	output :
	path '*.sexcheck'

	script:
	"""
	source /cm/shared/BCCHR-apps/env_vars/unset_BCM.sh
	source /cvmfs/soft.computecanada.ca/config/profile/bash.sh
	module load plink/1.9b_5.2-x86_64

	plink --vcf ${vcf} --split-x 2699520 154931044 --make-bed --vcf-half-call m --out ${vcf.simpleName}_plink_x_splitted
	plink --bfile ${vcf.simpleName}_plink_x_splitted --recode --impute-sex 0.5 0.8 --out ${vcf.simpleName}_plink_sex
	"""
}
