// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// Step merging the 2 vcf files (including the MT variants called against the reference genome and the shifted reference genome) for each individual
// bcftools norm remove the variants that are duplicated within the merge files (variants that were called against both references)

process MT_haplocheck {

	tag "${MT_FilterOut_sites_vcf.simpleName}"

//        publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/MT/Sample/", mode: 'copy'

	input :
	file MT_FilterOut_sites_vcf
	val assembly
	val batch
	val run

	output :
	path '*_haplocheck', emit : file

	script :
	"""
	sample_name=\$(echo ${MT_FilterOut_sites_vcf} | cut -d _ -f 1)
	if [ -a $params.outdir_ind/${assembly}/*/${run}/MT/Sample_vcf/\${sample_name}_MT_merged_filtered_trimmed_filtered_sites.vcf.gz ]; then
		touch ${MT_FilterOut_sites_vcf.simpleName}_haplocheck
	else
		/mnt/common/SILENT/Act3/haplocheck/./haplocheck --out ${MT_FilterOut_sites_vcf.simpleName}_haplocheck ${MT_FilterOut_sites_vcf}
	fi
	"""
}


