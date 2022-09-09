// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// Step merging the 2 vcf files (including the MT variants called against the reference genome and the shifted reference genome) for each individual
// bcftools norm remove the variants that are duplicated within the merge files (variants that were called against both references)

process MT_MergeVcfs {
        tag "${MT_call_variants.simpleName}"

        input :
        file MT_Liftover
	file MT_call_variants
	val assembly
	val batch
	val run

        output :
        path '*_MT_merged_uncollapsed.vcf.gz', emit : vcf

        script :
        """
        echo ${MT_call_variants.simpleName}
	sample_name=\$(echo ${MT_call_variants.simpleName} | cut -d _ -f 1)
	echo \$sample_name

	if [ -a $params.outdir_ind/${assembly}/*/${run}/MT/Sample/\${sample_name}_MT_merged_filtered_trimmed_filtered_sites.vcf.gz ]; then
		touch \${sample_name}_MT_merged_uncollapsed.vcf.gz
	else
        	gatk MergeVcfs \
		I=${MT_call_variants} \
		I=\${sample_name}_sorted_chrM_Homo_sapiens_assembly38_lifted_over.vcf \
		O=\${sample_name}_MT_merged_uncollapsed.vcf.gz
	fi
	"""
}
