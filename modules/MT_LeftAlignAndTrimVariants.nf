// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// MT. Filter merged Mutect2 calls


process MT_LeftAlignAndTrimVariants {
        tag "${MT_Filter_Mutect_Calls.simpleName}"

	input :
	file ref_genome_MT_file
	file ref_genome_MT_index
	file ref_genome_MT_dict
	file (MT_Filter_Mutect_Calls)
	file MT_Filter_Mutect_Calls_index
	val assembly
	val batch
	val run

	output :
	path '*_trimmed.vcf.gz', emit : vcf
	path '*_trimmed.vcf.gz.tbi', emit : index

	script :
	"""
	sample_name=\$(echo ${MT_Filter_Mutect_Calls} | cut -d _ -f 1)
	if [ -a $params.outdir_ind/${assembly}/*/${run}/MT/Sample_vcf/\${sample_name}_MT_merged_filtered_trimmed_filtered_sites.vcf.gz ]; then
		touch ${MT_Filter_Mutect_Calls.simpleName}_trimmed.vcf.gz
		touch ${MT_Filter_Mutect_Calls.simpleName}_trimmed.vcf.gz.tbi
	else
		gatk LeftAlignAndTrimVariants \
		-R Homo_sapiens_assembly38.chrM.fasta \
		-V ${MT_Filter_Mutect_Calls.simpleName}.vcf.gz \
		-O ${MT_Filter_Mutect_Calls.simpleName}_trimmed.vcf.gz \
		--split-multi-allelics \
		--dont-trim-alleles \
		--keep-original-ac
	fi
	"""
}

