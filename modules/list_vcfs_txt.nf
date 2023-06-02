// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// List the individuals files (vcf) that have been generated and that will be merged to obtain the aggregated dataset


process list_vcfs_txt {

        publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/${var_type}", mode: 'copy'
        
	input :
        file individual_vcf
	val assembly
	val batch
	val run
	val var_type 

        output :
        file '*.txt'

        script:
	if(var_type == "MT") {	
		"""
                find $params.outdir_ind/${assembly}/${batch}/${run}/${var_type}/Sample/ -name "*_filtered_sites.vcf.gz" > MT_vcfs.txt
		"""
	} else if (var_type == "SNV") {
                """
                find $params.outdir_ind/${assembly}/${batch}/${run}/${var_type}/Sample/ -name "*.g.vcf.gz" > ${var_type}_vcfs.txt
                """
	} else if (var_type == "SV") {
                """
                find $params.outdir_ind/${assembly}/${batch}/${run}/${var_type}/Sample/paragraph/ -name "*.vcf.gz" > ${var_type}_vcfs.txt
		"""
        } else {
        	"""
                find $params.outdir_ind/${assembly}/${batch}/${run}/${var_type}/Sample/ -name "*.vcf.gz" > ${var_type}_vcfs.txt
        	"""
	}

}

