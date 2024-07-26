// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developed to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// Quality check of the bam file with Picard CollectAlignmentSummaryMetrics
// This tool takes a SAM/BAM file input and produces metrics detailing the quality of the read alignments 
// as well as the proportion of the reads that passed machine signal-to-noise threshold quality filters. 

// June 2024 update by Kiana R: Nextflow's caching keeps track of already processed samples
// Thus no extra action is required to skip this process if the files have been already generated

process Picard_CollectAlignmentSummaryMetrics {
        tag "${bam}"

	publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/QC/Individuals/${bam.simpleName}/Picard_Metrics/", mode: 'copyNoFollow'

	input :
	file bam
	file bai 
	val assembly
	val batch
	val run
	
	output :
	file '*_Picard_Alignment*'

	script :
	"""
        gatk CollectAlignmentSummaryMetrics \
	--java-options "-Xmx2000M" \
	-I ${bam} \
	-O ${bam.simpleName}_Picard_Alignment
	"""
}

