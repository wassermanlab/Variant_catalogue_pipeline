// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// SNV calling. Call variants using deepvariant


process deepvariant_call {
        tag "${bam.simpleName}"

        publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/SNV/Sample/", mode: 'copyNoFollow'

	input :
	file reference
	file reference_index
	file bam
	file bai
	val assembly
	val batch
	val run

	output :
	path '*.g.vcf.gz', emit : deepvariant_gvcf
	path '*.vcf.gz', emit : deepvariant_vcf
	path '*.vcf.gz.tbi', emit : deepvariant_vcf_index

	script:
	"""
	if [ -a $params.outdir_ind/${assembly}/*/${run}/SNV/Sample/${bam.simpleName}.g.vcf.gz ]; then
		deepvariant_gvcf=\$(find $params.outdir_ind/${assembly}/*/${run}/SNV/Sample/ -name ${bam.simpleName}.g.vcf.gz) 
		deepvariant_vcf=\$(find $params.outdir_ind/${assembly}/*/${run}/SNV/Sample/ -name ${bam.simpleName}.vcf.gz)
		deepvariant_vcf_index=\$(find $params.outdir_ind/${assembly}/*/${run}/SNV/Sample/ -name ${bam.simpleName}.vcf.gz.tbi)
		ln -s \$deepvariant_gvcf .
		ln -s \$deepvariant_vcf .
		ln -s \$deepvariant_vcf_index .
	else
		/opt/deepvariant/bin/run_deepvariant \
		--num_shards=${task.cpus} \
		--intermediate_results_dir . \
		--model_type=WGS \
		--ref=${reference} \
		--reads=${bam} \
		--output_gvcf=${bam.simpleName}.g.vcf.gz \
		--output_vcf=${bam.simpleName}.vcf.gz
	fi
	"""
}
