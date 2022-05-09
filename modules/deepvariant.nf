// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// SNV calling. Call variants using deepvariant


process deepvariant_call {
        tag "${bam.simpleName}"

        publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/SNV/Sample/", mode: 'copy'

	input :
	file reference
	file reference_index
	file bam
	file bai
	val assembly
	val batch
	val run

	output :
	path '*_sorted.g.vcf.gz', emit : deepvariant_gvcf
	path '*_sorted.vcf.gz', emit : deepvariant_vcf
	path '*_sorted.vcf.gz.tbi', emit : deepvariant_vcf_index

	script:
	"""
	/opt/deepvariant/bin/run_deepvariant \
	--model_type=WGS \
	--ref=${reference} \
	--reads=${bam.simpleName}.bam \
	--output_gvcf=${bam.simpleName}.g.vcf.gz \
	--output_vcf=${bam.simpleName}.vcf.gz
	"""
}
