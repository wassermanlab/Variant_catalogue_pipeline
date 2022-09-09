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
	publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/MT/", mode: 'copy', pattern : 'MT.vcf.bgz*'
        publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/vcf_post_hail/", mode: 'copy', pattern : '*MT_filtered_with_geno*'
        publishDir "$params.outdir_pop/${assembly}/${run}/MT/Vcf_pre_annotation/", mode: 'copy', pattern : 'MT_filtered_frequ_only.vcf'


	input :
	file MT_Step1_input_tsv
	file MT_Step2_participant_data
	path MT_participants_to_subset
	file MT_Step3_participant_data
	val assembly
	val batch
	val run 
	file pon_prediction_table
	file artifact_prone_sites_bed
	file GRCh38_MT_local_fasta
	file GRCh38_MT_local_fai
	file mitotip_predictions_table

	output :
	path 'MT_post_hail_combined_sites_only.vcf.bgz', emit : vcf
	path 'MT_post_hail_combined_sites_only.vcf.bgz.tbi', emit : vcf_index
	path 'MT_filtered_frequ_only.vcf', emit : Hail_MT_frequ_only
	path 'sample_annotations_MT.txt', emit : sample_annotations
	path 'MT_stats_pass.txt', emit : stats_pass
	path 'MT_stats.txt', emit : stats

	script:
	"""
        #!/usr/bin/env python ../../../modules/Hail_variant_MT_QC.py $MT_Step1_input_tsv $MT_Step2_participant_data $MT_participants_to_subset $MT_Step3_participant_data $pon_prediction_table $artifact_prone_sites_bed $GRCh38_MT_local_fasta $GRCh38_MT_local_fai $mitotip_predictions_table
	"""
}


