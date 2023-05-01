// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// QC on the individual vcf using BcfTools stats
// For both the SNV and MT vcfs
// Another QC with the same tool is done on the aggregated vcf to do a different sample QC.

process Bcftools_stats {
	tag "${vcf}"

	publishDir "$params.outdir_pop/${assembly}/${run}/QC/${vcf.simpleName}/", mode: 'copy'

	input :
	file vcf
        file index
	val assembly
	val run

	output :
	file '*'

	script :
	"""
	ANNOTATEVARIANTS_INSTALL=/mnt/common/WASSERMAN_SOFTWARE/AnnotateVariants/
	source \$ANNOTATEVARIANTS_INSTALL/opt/miniconda3/etc/profile.d/conda.sh
	conda activate \$ANNOTATEVARIANTS_INSTALL/opt/AnnotateVariantsEnvironment

	bcftools stats ${vcf} > ${vcf.simpleName}_bcftools_stat
	"""
}
