// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developed to build the IBVL, a background variant library

// Overview of the process goal and characteristics:
// Alignment. fastq alignment with bwa mem 
// Sort and index with samtools

// Last edit: August 6, 2024 by Stephanie Petrone
// The bwa mem option -k 23 was added (changed from default of 19). 
// This was benchmarked with hap.py and had marginally better results; 
// Additionally, based on the HG002-4 samples, the alignment was faster. 


process align_sort_output_bam {
	label 'conda_annotate'
	tag "$sampleId"

	publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/BAM/"

	input :
	path reference
	path reference_index
	tuple val(sampleId), path ('read_pairs_ch')
	val assembly
	val batch
	val run

	output :
	path '*.bam', emit: samples_bam
	path '*.bam.bai', emit: samples_bam_index

	script:
	"""
	bwa mem -t ${task.cpus} -k 23 -R '@RG\\tID:${sampleId}\\tSM:${sampleId}' ${reference} ${read_pairs_ch} | samtools view -Sb | samtools sort -@ ${task.cpus} -o ${sampleId}_sorted.bam
	samtools index -@ ${task.cpus} ${sampleId}_sorted.bam
	"""
}

