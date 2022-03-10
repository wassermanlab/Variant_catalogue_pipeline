// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// BCFtools stats on the aggregated vcf (containing several samples)
// Part of the output (Part starting by PSC) is then loaded in R to obtain graphs with :
//	Number of singletons by individual
//	Number of variants by individuals
//	Het?hom ratio by individuals
// This values would allow to identify samples outside of the general distribution, which could be sign of low quality or contamination and should be removed from the final database

process count_bcftools_stats {
        tag "${vcf}"

	publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/QC/Out_MultiQC/variant_count/", mode: 'copy'

	input :
	file vcf 
	file vcf_index
	val assembly
	val batch
	val run
	
	output :
	file 'bcftools_stats.tsv'

	script :
	"""
	ANNOTATEVARIANTS_INSTALL=/mnt/common/WASSERMAN_SOFTWARE/AnnotateVariants/
        source \$ANNOTATEVARIANTS_INSTALL/opt/miniconda3/etc/profile.d/conda.sh
        conda activate \$ANNOTATEVARIANTS_INSTALL/opt/AnnotateVariantsEnvironment

        bcftools stats -s "-" ${vcf} > file.vchk 
	
	#Subset only the info wanted from the large table
	grep '^PSC' file.vchk > bcftools_stats.tsv
	"""
}

