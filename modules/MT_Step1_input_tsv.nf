// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// List the individuals files (vcf) that have been generated and that will be merged to obtain the aggregated dataset


process MT_Step1_input_tsv {

//        publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/${var_type}", mode: 'copy'
        
	input :
        file Sample_MT_Step1_input_tsv
	val assembly
	val batch
	val run

        output :
        file '*.tsv'

        script:
	"""
	cat $Sample_MT_Step1_input_tsv > MT_Step1_input_tsv.tsv
	"""
}

