process shift_back {
	tag "${MT_shifted_CollectMetrics}"

	publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/MT/QC/", mode: 'copyNoFollow', pattern: "*per_base_coverage.tsv"

	input :
	path MT_shifted_CollectMetrics
	path MT_CollectMetrics
	val assembly
	val batch
	val run
	
	output :
	path '*per_base_coverage.tsv', emit : per_base_coverage
	path '*_MT_Step1_input_tsv.tsv', emit : Sample_MT_Step1_input_tsv

	script:
	"""
        sample_name=\$(echo ${MT_shifted_CollectMetrics.simpleName} | sed 's/_.*//' )
	if [ -a $params.outdir_ind/${assembly}/${batch}/${run}/MT/QC/\${sample_name}_per_base_coverage.tsv]; then
		per_base_coverage=\$(find $params.outdir_ind/${assembly}/*/${run}/MT/QC/ -name \${sample_name}_per_base_coverage.tsv)
		ln -s \$per_base_coverage .
		touch \${sample_name}_MT_Step1_input_tsv.tsv
	else
		source /cm/shared/BCCHR-apps/env_vars/unset_BCM.sh
		source /cvmfs/soft.computecanada.ca/config/profile/bash.sh
		module load StdEnv/2020
		module load r/4.1.2

		Silent_Genomes_R=/mnt/common/SILENT/Act3/R/
		mkdir -p \${Silent_Genomes_R}/.local/R/\$EBVERSIONR/
		export R_LIBS=\${Silent_Genomes_R}/.local/R/\$EBVERSIONR/

		Rscript ../../../modules/shift_back.R $MT_shifted_CollectMetrics \${sample_name}_sorted_chrM_Homo_sapiens_assembly38_collect_wgs_metrics_non_control_region.chrM.interval_list.tsv
		mv per_base_coverage.tsv \${sample_name}_per_base_coverage.tsv	

		echo "\${sample_name}\t$params.outdir_ind/${assembly}/${batch}/${run}/MT/QC/\${sample_name}_per_base_coverage.tsv\t\${sample_name}" > \${sample_name}_MT_Step1_input_tsv.tsv
	fi
	"""
}
