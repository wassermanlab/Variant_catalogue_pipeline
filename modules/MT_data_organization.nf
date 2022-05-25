// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// Run a R script that calculates the distribution of Mitochondrial DNA Heteroplasmy for each variant in the cohort
// This is necessary to display the heteroplasmy histogram for each mitochondrial variant

process MT_data_organization {

        publishDir "$params.outdir_pop/${assembly}/${run}/Oracle_table/mt_ibvl_frequencies", mode: 'copy', pattern: "mt_ibvl_frequencies.tsv"
        publishDir "$params.outdir_pop/${assembly}/${run}/Oracle_table/mts/", mode: 'copy', pattern: "mts.tsv"
        publishDir "$params.outdir_pop/${assembly}/${run}/Oracle_table/variants/", mode: 'copy', pattern: "variants_MT.tsv"
        publishDir "$params.outdir_pop/${assembly}/${run}/Oracle_table/variants_transcripts/", mode: 'copy', pattern: "variants_transcripts_MT.tsv"
        publishDir "$params.outdir_pop/${assembly}/${run}/Oracle_table/variants_consequences/", mode: 'copy', pattern: "variants_consequences_MT.tsv"
        publishDir "$params.outdir_pop/${assembly}/${run}/Oracle_table/variant_annotations/", mode: 'copy', pattern: "variant_annotations_MT.tsv"
        publishDir "$params.outdir_pop/${assembly}/${run}/Oracle_table/mt_gnomad_frequencies/", mode: 'copy', pattern: "mt_gnomad_frequencies.tsv"
        publishDir "$params.outdir_pop/${assembly}/${run}/Oracle_table/genes/", mode: 'copy', pattern: "genes_MT.tsv"
        publishDir "$params.outdir_pop/${assembly}/${run}/Oracle_table/transcripts/", mode: 'copy', pattern: "transcripts_MT.tsv"

	input :
	path gnomad_MT_frequ
	path Hail_output
	path MT_annot
	val assembly
	val run

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

	Rscript ../../../modules/MT_data_organization.R $assembly $gnomad_MT_frequ $Hail_output $MT_annot $run
	"""
}
