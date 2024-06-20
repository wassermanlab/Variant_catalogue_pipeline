// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developed to build the IBVL, a background variant library

// Overview of the process goal and characteristics:
// Alignment. fastq alignment with bwa mem 
//	      Sort and index with samtools

// Process should be skipped if bam file already generated
// June 2024 update by Kiana R: Nextflow's caching keeps track of already processed samples
// Thus no extra action is required to skip this process if the files have been already generated

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
	bwa mem -t ${task.cpus} -R '@RG\\tID:${sampleId}\\tSM:${sampleId}' ${reference} ${read_pairs_ch} | samtools view -Sb | samtools sort -@ ${task.cpus} -o ${sampleId}_sorted.bam
	samtools index -@ ${task.cpus} ${sampleId}_sorted.bam
	"""
}

