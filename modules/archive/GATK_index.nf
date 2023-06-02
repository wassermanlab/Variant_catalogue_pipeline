// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// Merge all the samples (without joint calling) and create a vcf files with all the calls for each participant
// Index the vcf file

process GATK_index {

        publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/${var_type}/", mode: 'copy'

        input :
        file list_vcf
	val assembly
	val batch
	val run
	val var_type

        output:
	path '*.vcf.gz.tbi', emit : index

        script :
        """
        gatk --java-options "-Xmx4G" \
	IndexFeatureFile \
        -I ${var_type}_${run}.vcf.gz
        """
}
