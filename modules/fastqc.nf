// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// Pre-alignment QC : Fastqc

process fastqc {
        tag "$sample"
	
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
	if [ -a ${params.outdir_ind}/${assembly}/*/${run}/QC/Individuals/${sample}_sorted/Fastqc/${sample}.R1_fastqc.zip ]; then
		fastqc_R1=\$(find ${params.outdir_ind}/${assembly}/*/${run}/QC/Individuals/${sample}_sorted/Fastqc/ -name ${sample}.R1_fastqc.zip )
		fastqc_R2=\$(find ${params.outdir_ind}/${assembly}/*/${run}/QC/Individuals/${sample}_sorted/Fastqc/ -name ${sample}.R2_fastqc.zip )
		ln -s \$fastqc_R1 .
		ln -s \$fastqc_R2 .
	else
		fastqc -t ${task.cpus} ${reads.get(0)} ${reads.get(1)}
        fi
	"""
}
