// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// Run a python script with Hail to identify outliers variants
// It includes a step to filter the variants
// Future update : Include gnomAD frequency to each varaint, annotate the varaints using vep

process Hail_SV_QC {

        publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/QC/Aggregated/Hail/Variants/SV/", mode: 'copy', pattern : '*.html'
        publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/QC/Aggregated/Hail/Variants/SV/", mode: 'copy', pattern : 'SV_QC_report.txt'
        publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/vcf_post_hail/", mode: 'copy', pattern : 'SV_filtered_with_geno*'
        publishDir "$params.outdir_pop/${assembly}/${run}/SV/Vcf_pre_annotation/", mode: 'copy', pattern : 'SV_filtered_frequ_only*'

	input :
	file vcf_sample_filtered
	path sample_sex_file
	val assembly
	val batch
	val run 

	output :
	path '*.html', emit : graph
	path 'SV_QC_report.txt', emit : SV_QC_report

	path 'SV_filtered_with_geno.vcf.bgz*', emit : SV_filtered_with_geno
        path 'SV_filtered_frequ_only.vcf.bgz', emit : vcf_SV_filtered_frequ_only
        path 'SV_filtered_frequ_only.vcf.bgz.tbi', emit : index_SV_filtered_frequ_only

	conda '/home/BCRICWH.LAN/Solenne.Correard/miniconda3/envs/hail'

	script:
	"""
        #!/usr/bin/env python ../../../modules/Hail_SV_QC.py $vcf_sample_filtered $sample_sex_file
	"""
}


