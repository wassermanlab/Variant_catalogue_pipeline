// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics:
// Alignment. fastq alignment with bwa mem 
//	      Sort and index with samtools


process sorted_bam_files {
	tag "$sampleId"

	publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/BAM/", mode: 'copy'

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
	ANNOTATEVARIANTS_INSTALL=/mnt/common/WASSERMAN_SOFTWARE/AnnotateVariants/
	source \$ANNOTATEVARIANTS_INSTALL/opt/miniconda3/etc/profile.d/conda.sh
	conda activate \$ANNOTATEVARIANTS_INSTALL/opt/AnnotateVariantsEnvironment

	bwa mem -t 8 -R '@RG\\tID:${sampleId}\\tSM:${sampleId}' ${reference} ${read_pairs_ch} | samtools view -Sb | samtools sort -o ${sampleId}_sorted.bam	
	samtools index ${sampleId}_sorted.bam
	"""
}

