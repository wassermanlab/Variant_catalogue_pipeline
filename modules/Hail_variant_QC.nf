// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// Run a python script with Hail to identify outliers variants
// It includes a step to filter the variants
// Future update : Include gnomAD frequency to each varaint, annotate the varaints using vep

process Hail_variant_QC {

        publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/QC/Aggregated/Hail/Variants/", mode: 'copy', pattern : '*.html'
        publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/vcf_post_hail/", mode: 'copy', pattern : '*filtered_samples_variants.vcf.bgz*'
        publishDir "$params.outdir_pop/${assembly}/${run}/${var_type}/Vcf_pre_annotation/", mode: 'copy', pattern : 'SNV_filtered_frequ_*'

	input :
	file vcf_sample_filtered
	path sample_sex_file
	val assembly
	val batch
	val run 

	output :
	path '*.html', emit : graph

	path 'SNV_filtered_samples_variants*', emit : SNV_filtered_variants_ind_geno

	path 'SNV_filtered_frequ_total_xx.vcf.bgz', emit : SNV_filtered_variants_frequ_tot_xx
        path 'SNV_filtered_frequ_total_xx.vcf.bgz.tbi', emit : SNV_filtered_variants_frequ_tot_xx_index

	path 'SNV_filtered_tot_XX_XY.tsv.bgz', emit : SNV_frequ_tot_xx_xy_tsv



	conda '/home/BCRICWH.LAN/Solenne.Correard/miniconda3/envs/hail'

	script:
	"""
        #!/usr/bin/env python ../../../modules/Hail_variant_QC.py $vcf_sample_filtered $sample_sex_file
	"""
}


