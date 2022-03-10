// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// Run a R script to
//	Attribute a sex to each sample based on F coefficient, X and Y normalized coverage
//	Create graph with number of SNPs per individual, number of singletons per individual, het/hom ratio by individuals. 
// 		Blue lines in the graph are gnomAD threashold for sample to pass QC (Sample outside of the threasholds are removed in gnomAD)
//		No hard filter has been included to remove samples in the IBVL yet as the threasholds may differ from the threaholds applied by gnomAD
//		Visualization of the distributions using the graph will allow to define threasholds relevant to the IBVL


process sample_QC {

	publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/QC/Out_MultiQC/graph/", mode: 'copy'

	input :
	file plink_sex_inference
	val assembly
	val batch
	val run
	file singleton
	file SNPs_count
	file bcftools_stat
	file mosdepth 

	output :
	path '*.pdf', emit : graph
	path 'QC_sample.tsv', emit : sample_QC_file

	script:
	"""
	source /cm/shared/BCCHR-apps/env_vars/unset_BCM.sh
	source /cvmfs/soft.computecanada.ca/config/profile/bash.sh
	module load StdEnv/2020
	module load r/4.1.2

	Silent_Genomes_R=/mnt/common/SILENT/Act3/R/
	mkdir -p \${Silent_Genomes_R}/.local/R/\$EBVERSIONR/
	export R_LIBS=\${Silent_Genomes_R}/.local/R/\$EBVERSIONR/

	Rscript ../../../modules/sample_QC.R $assembly $batch $run $plink_sex_inference $singleton $SNPs_count $bcftools_stat
	"""
}
