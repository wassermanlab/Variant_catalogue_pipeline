// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// MT. Extract the MT reads from the bam file


process Extract_MT_Read {
        tag "${bam.simpleName}"

        input :
        file bam
	file bai
	val Mitochondrial_chromosome
	val assembly
	val batch
	val run


        output :
        file '*_chrM.bam'

        script :
        """
	sample_name=\$(echo ${bam.simpleName} | cut -d _ -f 1)
	if [ -a $params.outdir_ind/${assembly}/*/${run}/MT/Sample/\${sample_name}_MT_merged_filtered_trimmed_filtered_sites.vcf.gz ]; then
		touch \${sample_name}_chrM.bam
	else
		gatk PrintReads \
        	-L ${Mitochondrial_chromosome} \
        	--read-filter MateOnSameContigOrNoMappedMateReadFilter \
        	--read-filter MateUnmappedAndUnmappedReadFilter \
        	-I ${bam.simpleName}.bam \
        	--read-index ${bam.simpleName}.bam.bai \
        	-O ${bam.simpleName}_chrM.bam
        fi
	"""
}
