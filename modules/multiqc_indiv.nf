// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// MultiQC to aggregate all the individual QC data


process multiqc_indiv {

        publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/QC/Aggregated/multiqc/", mode: 'copy'

	input :
	file '*' 
	val assembly
	val batch
	val run

	output :
        file "multiqc_report.html"
        file "multiqc_data"

	script :
	"""
	multiqc $params.outdir_ind/${assembly}/${batch}/${run}/QC/
	"""
}
