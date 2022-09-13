// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// Run a R script that organize the SNV variants information in the tables expected to be displayed in the IBVL interface

process STR_data_organization {

        publishDir "$params.outdir_pop/${assembly}/${run}/Oracle_table/str/", mode: 'copy', pattern: "str.tsv"
	publishDir "$params.outdir_pop/${assembly}/${run}/Oracle_table/variants/", mode: 'copy', pattern: "variants_str.tsv"
        publishDir "$params.outdir_pop/${assembly}/${run}/Oracle_table/genes/", mode: 'copy', pattern: "genes_str.tsv"
        publishDir "$params.outdir_pop/${assembly}/${run}/Oracle_table/sv_consequences/", mode: 'copy', pattern: "sv_consequences_str.tsv"


	input :
	path STR_vcf
	path STR_catalogue
	val assembly
	val run
	val var_type 

	output :
	path '*'

	script:
	"""
	source /cm/shared/BCCHR-apps/env_vars/unset_BCM.sh
	source /cvmfs/soft.computecanada.ca/config/profile/bash.sh
	module load StdEnv/2020
	module load r/4.1.2

	Silent_Genomes_R=/mnt/common/SILENT/Act3/R/
	mkdir -p \${Silent_Genomes_R}/.local/R/\$EBVERSIONR/
	export R_LIBS=\${Silent_Genomes_R}/.local/R/\$EBVERSIONR/

	Rscript ../../../modules/STR_data_organization.R $assembly ${STR_vcf} ${STR_catalogue}  $run ${var_type}
	"""
}
