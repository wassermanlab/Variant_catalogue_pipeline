// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developed to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// Quality check of the bam file with Picard CollectWgsMetrics
// This tool collects metrics about the fractions of reads that pass base- and mapping-quality filters as well as coverage (read-depth) levels for WGS analyses.

// June 2024 update by Kiana R: Nextflow's caching keeps track of already processed samples
// Thus no extra action is required to skip this process if the files have been already generated

process Picard_CollectWgsMetrics {
        tag "${bam.simpleName}"
 
	publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/QC/Individuals/${bam.simpleName}/Picard_Metrics/", mode: 'copyNoFollow'

	input :
	file bam
	file bai
	path ref_genome 
	path reference_index
	val assembly
	val batch
	val run

	output :
	file '*_collect_wgs_metrics.txt' 

	script :
	"""
	gatk CollectWgsMetrics \
	--java-options "-Xmx8G" \
	-I ${bam} \
	-O ${bam.simpleName}_collect_wgs_metrics.txt \
	-R ${ref_genome}
	"""
}
