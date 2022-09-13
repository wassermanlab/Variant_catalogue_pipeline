// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// Run a R script that organize the SNV variants information in the tables expected to be displayed in the IBVL interface

process SNV_data_organization {
        tag "${SNV_annot_merged}"

	publishDir "$params.outdir_pop/${assembly}/${run}/Oracle_table/genomic_ibvl_frequencies/", mode: 'copy', pattern: "genomic_ibvl_frequencies_*"
        publishDir "$params.outdir_pop/${assembly}/${run}/Oracle_table/genomic_gnomad_frequencies/", mode: 'copy', pattern: "genomic_gnomad_frequencies_*"
        publishDir "$params.outdir_pop/${assembly}/${run}/Oracle_table/snvs/", mode: 'copy', pattern: "snvs_*"
        publishDir "$params.outdir_pop/${assembly}/${run}/Oracle_table/genes/", mode: 'copy', pattern: "genes_*"
        publishDir "$params.outdir_pop/${assembly}/${run}/Oracle_table/transcripts/", mode: 'copy', pattern: "transcripts_*"
        publishDir "$params.outdir_pop/${assembly}/${run}/Oracle_table/variants_transcripts/", mode: 'copy', pattern: "variants_transcripts_*"
        publishDir "$params.outdir_pop/${assembly}/${run}/Oracle_table/variants_consequences/", mode: 'copy', pattern: "variants_consequences_*"
        publishDir "$params.outdir_pop/${assembly}/${run}/Oracle_table/variants_annotations/", mode: 'copy', pattern: "variants_annotations_*"
        publishDir "$params.outdir_pop/${assembly}/${run}/Oracle_table/variants/", mode: 'copy', pattern: "variants_[A-Z0-9].tsv"
	publishDir "$params.outdir_pop/${assembly}/${run}/Oracle_table/variants/", mode: 'copy', pattern: "variants_[A-Z0-9][A-Z0-9].tsv"

	input :
	path gnomad_SNV_frequ
	path SNV_annot_merged
        path SNV_vcf
	val assembly
	val run
	path severity_table

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

	chr=\$(echo ${SNV_annot_merged.simpleName} | sed 's/^.*_\\([^_]*\\)\$/\\1/' )

	Rscript ../../../modules/SNV_data_organization.R $assembly gnomad_frequency_table_\${chr}.tsv ${SNV_vcf} ${SNV_annot_merged} $severity_table
	"""
}

//	vcf_name=\$(echo ${SNV_vcf.simpleName} | sed 's/_[^_]*\$//' )
//	chr=\$(echo ${SNV_vcf.simpleName} | sed 's/^.*_\\([^_]*\\)\$/\\1/' )


