// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// Run a R script to
//	Attribute a sex to each sample based on F coefficient, X and Y normalized coverage
//	Create graph with number of SNPs per individual, number of singletons per individual, het/hom ratio by individuals. 
// 		Blue lines in the graph are gnomAD threashold for sample to pass QC (Sample outside of the threasholds are removed in gnomAD)
//		No hard filter has been included to remove samples in the IBVL yet as the threasholds may differ from the threaholds applied by gnomAD
//		Visualization of the distributions using the graph will allow to define threasholds relevant to the IBVL


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
	path 'filtered_samples.vcf.bgz', emit : vcf_sample_filtered

	conda '/home/BCRICWH.LAN/Solenne.Correard/miniconda3/envs/hail'

	script:
	"""
        #!/usr/bin/env python ../../../modules/Hail_sample_QC.py $SNV_vcf
	"""
}
