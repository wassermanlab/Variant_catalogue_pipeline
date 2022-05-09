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
        publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/vcf_post_hail/", mode: 'copy', pattern : '*filtered_samples_variants.vcf.bgz'

	input :
	file vcf_sample_filtered
	val assembly
	val batch
	val run 

	output :
	path '*.html', emit : graph
	path 'filtered_samples_variants.vcf.bgz', emit : vcf_samples_variants_filtered

	conda '/home/BCRICWH.LAN/Solenne.Correard/miniconda3/envs/hail'

	script:
	"""
        #!/usr/bin/env python ../../../modules/Hail_variant_QC.py $vcf_sample_filtered
	"""
}
