// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// MT. Identify Duplicate Reads using MarkDuplicates
// This step identifies and tags duplicate reads in the aligned BAM files.



process MarkDuplicates {
	 tag "${bam_MT.baseName}"

        input :
        file bam_MT
	file bai_MT
	val assembly
	val batch
	val run

        output :
        path '*marked_duplicates.bam', emit : bam

        script:
        """
	sample_name=\$(echo ${bam_MT.baseName} | cut -d _ -f 1)
	if [ -a $params.outdir_ind/${assembly}/*/${run}/MT/Sample/\${sample_name}_MT_merged_filtered_trimmed_filtered_sites.vcf.gz ]; then
		touch ${bam_MT.baseName}_marked_duplicates.bam
	else
		gatk MarkDuplicates \
		I=${bam_MT.baseName}.bam \
		O=${bam_MT.baseName}_marked_duplicates.bam \
		M=${bam_MT.baseName}_marked_duplicates_metrics.txt
        fi
	"""
}
