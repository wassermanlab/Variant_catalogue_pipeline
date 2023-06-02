// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// Calculates the Transition / Transversion ratio as a function of SNP quality threshold. Only uses bi-allelic SNPs. The resulting output file has the suffix ".TsTv.qual".

process Vcftools_TsTv_by_qual {
	tag "${vcf}"

	publishDir "$params.outdir_pop/${assembly}/${run}/QC/${vcf.simpleName}/", mode: 'copy'

        input :
        file vcf
	file index
	val assembly
	val run

        output :
        file '*'

        script :
        """
        # Unload bcchr, and load cvmfs
        # unload_bcchr
        source /cm/shared/BCCHR-apps/env_vars/unset_BCM.sh
        # load cvmfs
        source /cvmfs/soft.computecanada.ca/config/profile/bash.sh

        module load StdEnv/2020
        module load vcftools

        vcftools --gzvcf ${vcf} --TsTv-by-qual --out ${vcf.simpleName}_Vcftools_TsTv_qual
	"""
}
