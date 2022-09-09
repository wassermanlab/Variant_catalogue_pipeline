// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// List the individuals files (vcf) that have been generated and that will be merged to obtain the aggregated dataset


process MT_Step3_metadata {
        
	input :
        file MT_Step3_metadata_sample
	val assembly
	val batch
	val run

        output :
        path '*.tsv'

        script:
	"""
	echo "entity:participant_id\ts\tcontamination\twgs_mean_coverage\tmt_mean_coverage" > header_MT_metadata_step3
	cat header_MT_metadata_step3 $MT_Step3_metadata_sample > MT_Step3_participant_data.tsv
	"""
}

