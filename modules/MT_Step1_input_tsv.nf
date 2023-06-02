// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// List the individuals files (vcf) that have been generated and that will be merged to obtain the aggregated dataset


process MT_Step1_input_tsv {
        
	input :
        file Sample_MT_Step1_input_tsv
	val assembly
	val batch
	val run

        output :
        file '*.tsv'

        script:
	"""
	sample_name=\$(echo ${Sample_MT_Step1_input_tsv.simpleName} | cut -d _ -f 1)
	if [ -a $params.outdir_ind/${assembly}/*/${run}/MT/Sample/\${sample_name}_MT_merged_filtered_trimmed_filtered_sites.vcf.gz ]; then
		touch MT_Step1_input_tsv.tsv
	else
		cat $Sample_MT_Step1_input_tsv > MT_Step1_input_tsv.tsv
	fi
	"""
}

