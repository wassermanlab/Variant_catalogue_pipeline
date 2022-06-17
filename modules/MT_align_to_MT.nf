// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// MT. Align the reads against the ref genome and the shifted ref genome


process align_to_MT {
        tag "$fastqfromsam.simpleName"

        input :
        path ref_genome_MT_file
	path ref_genome_MT_file_index
	file fastqfromsam
	val assembly
	val batch
	val run

        output :
        path '*.bam', emit : align_to_MT_bam
        path '*.bai', emit : align_to_MT_bai

        script:
	"""
	sample_name=\$(echo ${fastqfromsam.simpleName} | cut -d _ -f 1)
	if [ -a $params.outdir_ind/${assembly}/*/${run}/MT/Sample_vcf/\${sample_name}_MT_merged_filtered_trimmed_filtered_sites.vcf.gz ]; then
		touch ${fastqfromsam.baseName}_${ref_genome_MT_file.baseName}.bam
		touch ${fastqfromsam.baseName}_${ref_genome_MT_file.baseName}.bam.bai
	else
		ANNOTATEVARIANTS_INSTALL=/mnt/common/WASSERMAN_SOFTWARE/AnnotateVariants/
		source \$ANNOTATEVARIANTS_INSTALL/opt/miniconda3/etc/profile.d/conda.sh
		conda activate \$ANNOTATEVARIANTS_INSTALL/opt/AnnotateVariantsEnvironment

		bwa mem -R "@RG\\tID:${fastqfromsam.simpleName}\\tSM:${fastqfromsam.simpleName}\\tPL:illumina" ${ref_genome_MT_file} ${fastqfromsam.baseName}.fastq | samtools view -u -bS | samtools sort > ${fastqfromsam.simpleName}_${ref_genome_MT_file.baseName}.bam

        	samtools index ${fastqfromsam.baseName}_${ref_genome_MT_file.baseName}.bam
        fi
	"""
}
