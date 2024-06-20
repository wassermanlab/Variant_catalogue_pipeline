// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developed to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// Post-alignment QC
// Quality check of the bam file with Mosdepth

// June 2024 update by Kiana R: Nextflow's caching keeps track of already processed samples
// Thus no extra action is required to skip this process if the files have been already generated

process Mosdepth {
        tag "$bam"

	publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/QC/Individuals/${bam.simpleName}/Mosdepth/", mode: 'copy'

	input :
	file bam
	file bai
	val assembly
	val batch
	val run

	output : 
	path '*', emit : all_files
	path '*.mosdepth.summary.txt', emit : summary_stat

	script :
	"""
	mosdepth ${bam.simpleName} ${bam}
	"""
}
