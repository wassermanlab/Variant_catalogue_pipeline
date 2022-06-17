process MT_Step3_metadata_sample {
//	tag "${MT_shifted_CollectMetrics}"

//	publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/MT/QC/", mode: 'copy', pattern: "*per_base_coverage.tsv"

	input :
	path mosdepth
	path haplocheck
	val assembly
	val batch
	val run
	
	output :
	path '*', emit : MT_Step3_metadata_sample

	script:
	"""
        sample_name=\$(echo ${haplocheck.simpleName} | sed 's/_.*//' )
	if [ -a $params.outdir_ind/${assembly}/*/${run}/MT/Sample_vcf/\${sample_name}_MT_merged_filtered_trimmed_filtered_sites.vcf.gz ]; then
		touch \${sample_name}_conta_cov.tsv
	else
		source /cm/shared/BCCHR-apps/env_vars/unset_BCM.sh
		source /cvmfs/soft.computecanada.ca/config/profile/bash.sh
		module load StdEnv/2020
		module load r/4.1.2

		Silent_Genomes_R=/mnt/common/SILENT/Act3/R/
		mkdir -p \${Silent_Genomes_R}/.local/R/\$EBVERSIONR/
		export R_LIBS=\${Silent_Genomes_R}/.local/R/\$EBVERSIONR/

		Rscript ../../../modules/MT_Step3_metadata_sample.R \${sample_name}_sorted.mosdepth.summary.txt ${haplocheck}
		mv conta_cov.tsv \${sample_name}_conta_cov.tsv
	fi
	"""
}
