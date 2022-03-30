// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
//List the individuals files that will be merged to obtained the population mitochondrial calls


process list_vcfs_txt {

        publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/${var_type}", mode: 'copy'
        
	input :
        file expensionHunter_output
	val assembly
	val batch
	val run
	val var_type 

        output :
        file '*.txt'

        script:
	if(var_type == "MT") {	
		"""
                find $params.outdir_ind/${assembly}/${batch}/${run}/${var_type}/ -name "*_filtered_sites.vcf.gz" > MT_vcfs.txt
        	"""
	} else if (var_type == "SV") {
                """
                find $params.outdir_ind/${assembly}/${batch}/${run}/${var_type}/paragraph/ -name "*.vcf.gz" > ${var_type}_vcfs.txt
		"""
        } else {
        	"""
                find $params.outdir_ind/${assembly}/${batch}/${run}/${var_type}/ -name "*.vcf.gz" > ${var_type}_vcfs.txt
        	"""
	}

}

