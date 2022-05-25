// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// List the individuals files (vcf) that have been generated and that will be merged to obtain the aggregated dataset


process MT_Step2_participant_data {

//        publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/${var_type}", mode: 'copy'
        
	input :
        file Sample_MT_Step2_participant_data
	file Sample_list
	val assembly
	val batch
	val run

        output :
        path '*.tsv', emit : MT_Step2_participant_data_tsv
	path '*.txt', emit : participants_to_subset_txt

        script:
	"""
	echo "entity:participant_id\ts\tvcf_path" > header_MT_Step2_participant_data
	cat header_MT_Step2_participant_data $Sample_MT_Step2_participant_data > MT_Step2_participant_data.tsv
	
	echo "participant" > header_participants_to_subset
	cat header_participants_to_subset $Sample_list > MT_participants_to_subset.txt
	"""
}

