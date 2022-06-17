// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// Step merging the 2 vcf files (including the MT variants called against the reference genome and the shifted reference genome) for each individual
// bcftools norm remove the variants that are duplicated within the merge files (variants that were called against both references)

process MT_MergeVcfs {
        tag "${MT_call_variants.simpleName}"

//        publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/MT/Sample/", mode: 'copy'

        input :
        file MT_Liftover
	file MT_call_variants
	val assembly
	val batch
	val run

        output :
        path '*_merged.vcf.gz', emit : vcf
	path '*_merged.vcf.gz.tbi', emit : index

        script :
        """
        echo ${MT_call_variants.simpleName}
	sample_name=\$(echo ${MT_call_variants.simpleName} | cut -d _ -f 1)
	echo \$sample_name

	if [ -a $params.outdir_ind/${assembly}/*/${run}/MT/Sample_vcf/\${sample_name}_MT_merged_filtered_trimmed_filtered_sites.vcf.gz ]; then
		touch \${sample_name}_MT_merged.vcf.gz
		touch \${sample_name}_MT_merged.vcf.gz.tbi
	else
		singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/gatk4-4.2.0.sif \
        	gatk MergeVcfs \
		I=${MT_call_variants} \
		I=\${sample_name}_sorted_chrM_Homo_sapiens_assembly38_lifted_over.vcf \
		O=\${sample_name}_MT_merged_uncollapsed.vcf.gz

		ANNOTATEVARIANTS_INSTALL=/mnt/common/WASSERMAN_SOFTWARE/AnnotateVariants/
		source \$ANNOTATEVARIANTS_INSTALL/opt/miniconda3/etc/profile.d/conda.sh
		conda activate \$ANNOTATEVARIANTS_INSTALL/opt/AnnotateVariantsEnvironment

		bcftools norm --rm-dup both \${sample_name}_MT_merged_uncollapsed.vcf.gz -O z -o \${sample_name}_MT_merged.vcf.gz 
		bcftools index -t \${sample_name}_MT_merged.vcf.gz 
	fi
	"""
}
