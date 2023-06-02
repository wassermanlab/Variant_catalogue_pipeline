// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
//List the individuals files that will be merged to obtained the population mitochondrial calls


process MT_vcfs_txt {

        publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/Mutect2/", mode: 'copy'
        
	input :
        file MT_FilterOut
	val assembly
	val batch
	val run

        output :
        file '*.txt'

        script:
        """
        find $params.outdir_ind/${assembly}/${batch}/${run}/Mutect2/ -name "*_filtered_sites.vcf.gz" > MT_vcfs.txt
        """
}
