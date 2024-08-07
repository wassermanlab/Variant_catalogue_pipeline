// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developed to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// Pre-alignment QC : Fastqc


process fastqc {
        tag "$sample"
        container = "https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0"	
	publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/QC/Individuals/${sample}_sorted/Fastqc/", mode: 'copyNoFollow'

	input:
        tuple (val(sample), file(reads)) 
        path outdir_ind
	val assembly
	val batch
	val run

	output :
	file ("*.zip")

	script:
        """
	fastqc -t ${task.cpus} ${reads.get(0)} ${reads.get(1)}
	"""
}
