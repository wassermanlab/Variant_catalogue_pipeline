// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// MT. Call the variants against the ref genome and the shifted ref genome


process MT_call_variants {
        tag "${MarkDuplicates_bam_MT.baseName}"

	input :
	path ref_genome_MT_file
        path ref_genome_MT_file_index	
	path ref_genome_MT_file_dict
	file MarkDuplicates_bam_MT
        file MarkDuplicates_bam_MT_bai
	val Mitochondrial_chromosome
	val assembly
	val batch
	val run

	output :
	path '*_Mutect2.vcf.gz', emit: Mutect2_vcf
        path '*_Mutect2.vcf.gz.tbi', emit: Mutect2_vcf_index
	path '*.vcf.gz.stats', emit: Mutect2_stat

	script:
	"""
	sample_name=\$(echo ${MarkDuplicates_bam_MT} | cut -d _ -f 1)
	if [ -a $params.outdir_ind/${assembly}/*/${run}/MT/Sample/\${sample_name}_MT_merged_filtered_trimmed_filtered_sites.vcf.gz ]; then
		touch ${MarkDuplicates_bam_MT.baseName}_Mutect2.vcf.gz
		touch ${MarkDuplicates_bam_MT.baseName}_Mutect2.vcf.gz.tbi
		touch ${MarkDuplicates_bam_MT.baseName}_Mutect2.vcf.gz.stats
	else
		gatk Mutect2 \
		-R ${ref_genome_MT_file} \
		-I ${MarkDuplicates_bam_MT.baseName}.bam \
		-L chrM \
		--mitochondria-mode \
		--annotation StrandBiasBySample \
		--max-reads-per-alignment-start 75 \
		--max-mnp-distance 0 \
		-O ${MarkDuplicates_bam_MT.baseName}_Mutect2.vcf.gz
	fi
	"""
}
