// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// MT. Filter merged Mutect2 calls


process MT_Filter_Mutect_Calls {
        tag "${MT_MergeVcfs.simpleName}"

	input :
	file ref_genome_MT
	file ref_genome_MT_index
	file ref_genome_MT_dict
	file MT_MergeVcfs
	file MT_MergeVcfs_index
	file MT_MergeVcfs_stat
	val assembly
	val batch
	val run

	output :
	path '*_filtered.vcf.gz', emit : vcf
	path '*_filtered.vcf.gz.tbi', emit : index

	script :
	"""
	sample_name=\$(echo ${MT_MergeVcfs} | cut -d _ -f 1)
	if [ -a $params.outdir_ind/${assembly}/*/${run}/MT/Sample_vcf/\${sample_name}_MT_merged_filtered_trimmed_filtered_sites.vcf.gz ]; then
		touch ${MT_MergeVcfs.simpleName}_filtered.vcf.gz
		touch ${MT_MergeVcfs.simpleName}_filtered.vcf.gz.tbi
	else
        	gatk FilterMutectCalls \
		-V ${MT_MergeVcfs.simpleName}.vcf.gz \
		-R Homo_sapiens_assembly38.chrM.fasta \
		--stats ${MT_MergeVcfs.simpleName}.stats \
		--max-alt-allele-count 4 \
		--mitochondria-mode \
		-O ${MT_MergeVcfs.simpleName}_filtered.vcf.gz
	fi
	"""
}

