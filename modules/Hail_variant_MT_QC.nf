// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// Run a python script with Hail to identify outliers variants
// It includes a step to filter the variants
// Future update : Include gnomAD frequency to each varaint, annotate the varaints using vep

process Hail_variant_MT_QC {

        publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/QC/Aggregated/Hail/Variants/MT/", mode: 'copy', pattern : 'sample_annotations.txt'
        publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/QC/Aggregated/Hail/Variants/MT/", mode: 'copy', pattern : 'MT_stats.txt'
        publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/vcf_post_hail/", mode: 'copy', pattern : '*combined_sites_only.vcf.bgz*'

	input :
	file MT_Step1_input_tsv
	file MT_Step2_participant_data
	path MT_participants_to_subset
	file MT_Step3_participant_data
	val assembly
	val batch
	val run 

	output :
	path 'MT_post_hail_combined_sites_only.vcf.bgz', emit : vcf
	path 'MT_post_hail_combined_sites_only.vcf.bgz.tbi', emit : vcf_index
	path 'MT_Step3_output_dir/reduced_annotations.txt', emit : Hail_reduced_annotations
	path 'sample_annotations_MT.txt', emit : sample_annotations
	path 'MT_stats_pass.txt', emit : stats_pass
	path 'MT_stats.txt', emit : stats

	conda '/home/BCRICWH.LAN/Solenne.Correard/miniconda3/envs/hail'

	script:
	"""
        #!/usr/bin/env python ../../../modules/Hail_variant_MT_QC.py $MT_Step1_input_tsv $MT_Step2_participant_data $MT_participants_to_subset $MT_Step3_participant_data
	"""
}


