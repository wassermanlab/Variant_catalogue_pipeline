// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// Run a R script that organize the SNV variants information in the tables expected to be displayed in the IBVL interface

process MEI_data_organization {
        tag "${SV_vcf}"

	publishDir "$params.outdir_pop/${assembly}/${run}/Oracle_table/genomic_ibvl_frequencies/", mode: 'copy', pattern: "genomic_ibvl_frequencies_*"
        publishDir "$params.outdir_pop/${assembly}/${run}/Oracle_table/svs/", mode: 'copy', pattern: "svs_*"
        publishDir "$params.outdir_pop/${assembly}/${run}/Oracle_table/genes/", mode: 'copy', pattern: "genes_sv_*"
	publishDir "$params.outdir_pop/${assembly}/${run}/Oracle_table/variants/", mode: 'copy', pattern: "variants_sv_MEI_[A-Z0-9][A-Z0-9].tsv"
	publishDir "$params.outdir_pop/${assembly}/${run}/Oracle_table/sv_consequences/", mode: 'copy', pattern: "sv_consequences_*"

	input :
	path MEI_vcf
	path MEI_annot_merged
	val assembly
	val run
	val var_type 
	path sex_table

	output :
	path '*'

	script:
	"""
	source /cm/shared/BCCHR-apps/env_vars/unset_BCM.sh
	source /cvmfs/soft.computecanada.ca/config/profile/bash.sh
	module load StdEnv/2020
	module load r/4.1.2

	Silent_Genomes_R=/mnt/common/SILENT/Act3/R/
	mkdir -p \${Silent_Genomes_R}/.local/R/\$EBVERSIONR/
	export R_LIBS=\${Silent_Genomes_R}/.local/R/\$EBVERSIONR/

	vcf_name=\$(echo ${MEI_vcf.simpleName} | sed 's/_[^_]*\$//' )
	chr=\$(echo ${MEI_vcf.simpleName} | sed 's/^.*_\\([^_]*\\)\$/\\1/' )

	Rscript ../../../modules/MEI_data_organization.R $assembly ${MEI_vcf} \${vcf_name}_\${chr}_${var_type}_annotation_table_merged_nohash.tsv $sex_table $run ${var_type}
	"""
}
