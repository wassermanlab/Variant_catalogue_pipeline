// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// Run a python script with Hail to identify outliers variants
// It includes a step to filter the variants
// Future update : Include gnomAD frequency to each varaint, annotate the varaints using vep

process Hail_variant_QC {
        //  maxForks 1
        publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/QC/Aggregated/Hail/Variants/SNV/", mode: 'copy', pattern : '*.html'
        publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/QC/Aggregated/Hail/Variants/SNV/", mode: 'copy', pattern : 'SNV_indel_QC_report.txt'
        publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/vcf_post_hail/", mode: 'copy', pattern : 'SNV_filtered_with_geno*'
        publishDir "$params.outdir_pop/${assembly}/${run}/SNV/Vcf_pre_annotation/", mode: 'copy', pattern : 'SNV_filtered_frequ_only*'
		publishDir "$params.outdir_pop/${assembly}/${run}/SNV/Vcf_pre_annotation/large_indels/", mode: 'copy', pattern : 'SNV_large_indels*'
		publishDir "$params.outdir_pop/${assembly}/${run}/SNV/Vcf_pre_annotation/autosomal/", mode: 'copy', pattern : 'SNV_autosomal.vcf*'

	input :
	file vcf_sample_filtered
	path sample_sex_file
	val assembly
	val ref
	val ref_index
	val batch
	val run
	each var_qc_interval

	output :
	path '*.html', emit : graph
	path 'SNV_variant_QC_report*.txt', emit : SNV_QC_report

	path 'SNV_filtered_with_geno*.vcf.bgz*', emit : SNV_filtered_with_geno
    path 'SNV_filtered_frequ_only*.vcf.bgz', emit : vcf_SNV_filtered_frequ_only
    path 'SNV_filtered_frequ_only*.vcf.bgz.tbi', emit : index_SNV_filtered_frequ_only

	path 'SNV_large_indels*.vcf*'
	path '*.vcf*'

	script:
	"""
        #!/usr/bin/env python ../../../modules/Hail_variant_QC.py $vcf_sample_filtered $sample_sex_file $params.tmp_dir $assembly $ref $ref_index $var_qc_interval
	"""
}


