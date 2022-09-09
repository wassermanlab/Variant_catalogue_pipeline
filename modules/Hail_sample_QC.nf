// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// Run a python script involving Hail to create graph to identify outliers samples
// Includes a final filtering step 

process Hail_sample_QC {

	publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/QC/Aggregated/Hail/Samples/", mode: 'copy', pattern : '*.html'
	publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/vcf_post_hail/", mode: 'copy', pattern : '*filtered_samples.vcf.bgz'

	input :
	file SNV_vcf
	val assembly
	val batch
	val run 

	output :
	path '*.html', emit : graph
	path '*filtered_samples.vcf.bgz', emit : vcf_sample_filtered
	path '*filtered_samples_sex.tsv', emit : filtered_sample_sex

	script:
	"""
        #!/usr/bin/env python ../../../modules/Hail_sample_QC.py $SNV_vcf
	"""
}
