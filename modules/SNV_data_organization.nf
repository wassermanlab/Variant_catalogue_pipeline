// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// Run a R script that organize the SNV variants information in the tables expected to be displayed in the IBVL interface

process SNV_data_organization {
        tag "${SNV_vcf}"

        publishDir "$params.outdir_pop/${assembly}/${run}/Oracle_table/", mode: 'copy'

	input :
	path gnomad_SNV_frequ
	path SNV_vcf
	path SNV_annot_merged
	val assembly
	val run
	path sex_table

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

	vcf_name=\$(echo ${SNV_vcf.simpleName} | sed 's/_[^_]*\$//' )
	chr=\$(echo ${SNV_vcf.simpleName} | sed 's/^.*_\\([^_]*\\)\$/\\1/' )

	Rscript ../../../modules/SNV_data_organization.R $assembly gnomad_frequency_table_\${chr}.tsv ${SNV_vcf} \${vcf_name}_\${chr}_annotation_table_merged_nohash.tsv $sex_table $run
	"""
}
