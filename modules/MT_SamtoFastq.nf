// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// MT. Convert the bam file containing only the reads mapping against the MT chromosome into fastq

process MT_SamtoFastq {
        tag "${Extract_MT_Read.baseName}"

        input :
        file Extract_MT_Read

        output :
        path '*.fastq', emit : fastq_MT
//	path '*.fastq.fai', emit : fastq_MT_index

        script :
        """
        gatk SamToFastq \
        INPUT=${Extract_MT_Read.baseName}.bam \
        FASTQ=${Extract_MT_Read.baseName}.fastq \
        INTERLEAVE=true \
        NON_PF=true
	"""
}

