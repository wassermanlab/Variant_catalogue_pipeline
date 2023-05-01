// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// Run a python script with Hail to identify outliers variants
// It includes a step to filter the variants
// Future update : Include gnomAD frequency to each varaint, annotate the varaints using vep

process Hail_variant_QC {

        publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/QC/Aggregated/Hail/Variants/SNV/", mode: 'copy', pattern : '*.html'
        publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/QC/Aggregated/Hail/Variants/SNV/", mode: 'copy', pattern : 'SNV_indel_QC_report.txt'
        publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/vcf_post_hail/", mode: 'copy', pattern : 'SNV_filtered_with_geno*'
        publishDir "$params.outdir_pop/${assembly}/${run}/SNV/Vcf_pre_annotation/", mode: 'copy', pattern : 'SNV_filtered_frequ_only*'

	input :
	file vcf_sample_filtered
	path sample_sex_file
	val assembly
	val batch
	val run 

	output :
	path '*.html', emit : graph
	path 'SNV_indel_QC_report.txt', emit : SNV_QC_report

	path 'SNV_filtered_with_geno.vcf.bgz*', emit : SNV_filtered_with_geno
        path 'SNV_filtered_frequ_only.vcf.bgz', emit : vcf_SNV_filtered_frequ_only
        path 'SNV_filtered_frequ_only.vcf.bgz.tbi', emit : index_SNV_filtered_frequ_only

	script:
	"""
        #!/usr/bin/env python ../../../modules/Hail_variant_QC.py $vcf_sample_filtered $sample_sex_file $params.tmp_dir $assembly
	"""
}


