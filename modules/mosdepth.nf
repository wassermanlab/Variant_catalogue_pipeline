// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// Post-alignment QC
// Quality check of the bam file with Mosdepth

process Mosdepth {
        tag "$bam"

	publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/QC/Individuals/${bam.simpleName}/Mosdepth/", mode: 'copyNoFollow'

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
	if [ -a $params.outdir_ind/${assembly}/*/${run}/QC/Individuals/${bam.simpleName}/Mosdepth/${bam.simpleName}.mosdepth.global.dist.txt  ]; then
		glob_dist=\$(find $params.outdir_ind/${assembly}/*/${run}/QC/Individuals/${bam.simpleName}/Mosdepth/ -name ${bam.simpleName}.mosdepth.global.dist.txt)
		summary=\$(find $params.outdir_ind/${assembly}/*/${run}/QC/Individuals/${bam.simpleName}/Mosdepth/ -name ${bam.simpleName}.mosdepth.summary.txt)
		per_base_bed=\$(find $params.outdir_ind/${assembly}/*/${run}/QC/Individuals/${bam.simpleName}/Mosdepth/ -name ${bam.simpleName}.per-base.bed.gz)
		per_base_index=\$(find $params.outdir_ind/${assembly}/*/${run}/QC/Individuals/${bam.simpleName}/Mosdepth/ -name ${bam.simpleName}.per-base.bed.gz.csi)
		ln -s \$glob_dist .
		ln -s \$summary .
		ln -s \$per_base_bed  .
		ln -s \$per_base_index .
	else
		mosdepth ${bam.simpleName} ${bam}
	fi
	"""
}
