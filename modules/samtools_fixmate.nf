// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics:
// Alignment. fastq alignment with bwa mem 
//	      Sort and index with samtools


process samtools_fixmate {
	tag "${bam.SimpleName}"

//	publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/BAM/", mode: 'copy'

	input :
	file bam
	file bai
	val assembly
	val batch
	val run

	output :
	path '*_fixmate_ordered.bam', emit: samples_fixmate_bam
	path '*_fixmate_ordered.bam.bai', emit: samples_fixmate_bam_index

	script:
	"""
	ANNOTATEVARIANTS_INSTALL=/mnt/common/WASSERMAN_SOFTWARE/AnnotateVariants/
	source \$ANNOTATEVARIANTS_INSTALL/opt/miniconda3/etc/profile.d/conda.sh
	conda activate \$ANNOTATEVARIANTS_INSTALL/opt/AnnotateVariantsEnvironment

	# Resort the bam file by query name for samtools fixmate (coordiante-sorted bam does not work)
	samtools sort -n -O BAM -@ 20  ${bam} > ${bam.SimpleName}_nsorted.bam

	# Samtools fixmate will add MQ tags
	samtools fixmate -m -O BAM -@ 20  ${bam.SimpleName}_nsorted.bam  ${bam.SimpleName}_fixmate.bam

	# Now sort bam file by coordinates to resume the pipeline 
	samtools sort  -@ 20  ${bam.SimpleName}_fixmate.bam -o ${bam.SimpleName}_fixmate_ordered.bam

	# index the sorted bam file
	samtools index  -@ 20  ${bam.SimpleName}_fixmate_ordered.bam
	"""
}

