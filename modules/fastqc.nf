// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// Pre-alignment QC : Fastqc

process fastqc {
        tag "$sample"
	
	publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/QC/${sample}_sorted/Fastqc/", mode: 'copy'

	input:
        tuple (val(sample), file(reads)) 
        path outdir_ind
	val assembly
	val batch
	val run

	output :
//	file ("${sample}_fastqc/*.zip")
	file ("*.zip")

	script:
//        mkdir -p ${sample}_fastqc
//        fastqc --outdir ${sample}_fastqc -t ${task.cpus} ${reads.get(0)} ${reads.get(1)}
        """
	fastqc -t ${task.cpus} ${reads.get(0)} ${reads.get(1)}
        """
}
