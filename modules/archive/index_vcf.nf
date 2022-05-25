// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// Index compressed vcf

process index_vcf {
        publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/vcf_post_hail/", mode: 'copy'

	input :
	file vcf_file
	val assembly
	val batch
	val run

	output :
	path '*.vcf.bgz.tbi', emit : index
	

	script :
	"""
        singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/gatk4-4.2.0.sif \
	gatk --java-options "-Xmx4G" \
	IndexFeatureFile \
        -I ${vcf_file.simpleName}.vcf.bgz	
	"""
}
