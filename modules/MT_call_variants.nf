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

	output :
	path '*_Mutect2.vcf.gz', emit: Mutect2_vcf
        path '*_Mutect2.vcf.gz.tbi', emit: Mutect2_vcf_index
	path '*.vcf.gz.stats', emit: Mutect2_stat

	script:
	"""
	gatk Mutect2 \
	-R ${ref_genome_MT_file} \
	-I ${MarkDuplicates_bam_MT.baseName}.bam \
	-L chrM \
	--mitochondria-mode \
	--annotation StrandBiasBySample \
	--max-reads-per-alignment-start 75 \
	--max-mnp-distance 0 \
	-O ${MarkDuplicates_bam_MT.baseName}_Mutect2.vcf.gz
	"""
}
