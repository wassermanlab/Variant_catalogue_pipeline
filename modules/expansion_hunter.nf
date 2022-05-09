// Nextflow process
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// SV calling. Genotype the Short Tendem Repeat (STR) using Expension Hunter
// Rename the STR, compress the vcf and index the compressed vcf

// Future possible improvement : Use ExpensionHunterDeNovo


process expansion_hunter {
	tag "${bam.simpleName}"

        publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/STR/Sample/", mode: 'copy'

	input:
	file bam
	file bai
	file reference
	file reference_index
	file variant_catalog
        val assembly
        val batch
        val run

	output:
	path '*_str.vcf.gz', emit : vcf
        path '*_str.vcf.gz.tbi', emit : vcf_index

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
