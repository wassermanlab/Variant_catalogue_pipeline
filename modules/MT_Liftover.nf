// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// This step returns the variant calls back to the standard numbering system with the original alignment (OA) tags.

process MT_Liftover {
        tag "${MT_call_variants_shifted.simpleName}"

        publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/MT/QC/${MT_call_variants_shifted.simpleName}/Liftover/",  pattern: "*_rejected_variants.vcf",  mode: 'copyNoFollow'

	input :
        file MT_call_variants_shifted
        file MT_call_variants_shifted_index
	file ref_genome_MT_file
	file ref_genome_MT_file_dict
	file bwa_index_ref_genome
	file ShiftBack_chain_MT_file
	val assembly
	val batch
	val run

        output :
        path '*_lifted_over.vcf', emit : lifted_vcf
	path '*_rejected_variants.vcf', emit : rejected_vcf


        script :
        """
	sample_name=\$(echo ${MT_call_variants_shifted} | cut -d _ -f 1)
	if [ -a $params.outdir_ind/${assembly}/*/${run}/MT/Sample_vcf/\${sample_name}_MT_merged_filtered_trimmed_filtered_sites.vcf.gz ]; then
		touch ${MT_call_variants_shifted.simpleName}_lifted_over.vcf
		rejected_vcf=\$(find $params.outdir_ind/${assembly}/*/${run}/MT/QC/\${sample_name}*/Liftover/ -name \${sample_name}*_rejected_variants.vcf)
		ln -s \$rejected_vcf .
	else
		gatk LiftoverVcf \
		I=${MT_call_variants_shifted} \
		O=${MT_call_variants_shifted.simpleName}_lifted_over.vcf \
		CHAIN=${ShiftBack_chain_MT_file} \
		REJECT=${MT_call_variants_shifted.simpleName}_rejected_variants.vcf \
		R=${ref_genome_MT_file}
	fi
	"""
}
