// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// Count the number of singletons in each vcf.
// The number of singleton is used to filter samples based on quality.
// BCFtools_stats is also used to count the number of singletons but they give rather different numbers. A comparison is needed.
// The output is loaded into QC_Sample, the R script creatin the graph to determine the QC threasholds.

process count_variants_vcftools {
        tag "${vcf}"

	publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/QC/Aggregated/variant_count/", mode: 'copy'

	input :
	file vcf 
	file vcf_index
	val assembly
	val batch
	val run
	
	output :
	file 'singleton_per_ind_vcftools.tsv'

	script :
	"""
	# Unload bcchr, and load cvmfs
        # unload_bcchr
        source /cm/shared/BCCHR-apps/env_vars/unset_BCM.sh
        # load cvmfs
        source /cvmfs/soft.computecanada.ca/config/profile/bash.sh

        module load StdEnv/2020
        module load vcftools
	
	vcftools --singletons --gzvcf ${vcf}
	
	# Number of lines with each sample = number of singleton for each sample
        awk '{print \$5}' out.singletons | sort | uniq -c > singleton_per_ind_vcftools.tsv
	"""
}

