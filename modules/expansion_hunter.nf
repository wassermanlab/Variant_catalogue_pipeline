process expansion_hunter {
	tag "${bam.simpleName}"

//	publishDir "PipeRes", mode: 'copy'
        publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/STR", mode: 'copy'

	input:
//        tuple(val(sample_name), path(bam), path(index))
	file bam
	file bai
	file reference
	file reference_index
	file variant_catalog
        val assembly
        val batch
        val run

	output:
//	path '*.json', emit : json
	path '*_str.vcf.gz', emit : vcf
        path '*_str.vcf.gz.tbi', emit : vcf_index
//        tuple(file("${sample_name}_str.vcf"), val(sample_name))

	script:
	"""
	${params.ExpansionHunter_dir}/ExpansionHunter \
	--output-prefix ${bam.simpleName} \
	--reference $reference  \
	--reads ${bam} \
	--variant-catalog ${variant_catalog} \
	-n ${task.cpus} 
	
	# Unload bcchr, and load cvmfs
        # unload_bcchr
        source /cm/shared/BCCHR-apps/env_vars/unset_BCM.sh
        # load cvmfs
        source /cvmfs/soft.computecanada.ca/config/profile/bash.sh

        module load StdEnv/2020
	module load bcftools
		
	bcftools view -O z -o ${bam.simpleName}_str_noID.vcf.gz ${bam.simpleName}.vcf
	bcftools index ${bam.simpleName}_str_noID.vcf.gz
	bcftools annotate --set-id '%CHROM\\_%POS\\_%END\\_%REF\\_%ALT' -O z -o ${bam.simpleName}_str.vcf.gz ${bam.simpleName}_str_noID.vcf.gz
	bcftools index --tbi ${bam.simpleName}_str.vcf.gz
	"""
}
