// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// MultiQC to aggregate all the population QC data


process multiqc_pop {

	publishDir "$params.outdir_pop/${assembly}/${run}/QC/MultiQC_pop/${variant_type}/", mode: 'copy'

        input :
        file '*'
	val assembly
	val run
	val variant_type

        output :
        file '*'

        script :
        """
	multiqc $params.outdir_pop/${assembly}/${run}/QC/
        """
}

