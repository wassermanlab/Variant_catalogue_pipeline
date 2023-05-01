process shift_back {
	tag "${MT_shifted_CollectMetrics}"

	publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/MT/QC/", mode: 'copyNoFollow', pattern: "*per_base_coverage.tsv"

	input :
	path MT_shifted_CollectMetrics
	path MT_CollectMetrics
	val assembly
	val batch
	val run
	val path_R_libraries
	
	output :
	path '*per_base_coverage.tsv', emit : per_base_coverage
	path '*_MT_Step1_input_tsv.tsv', emit : Sample_MT_Step1_input_tsv

	script:
	"""
        sample_name=\$(echo ${MT_shifted_CollectMetrics.simpleName} | sed 's/_.*//' )
	if [ -a $params.outdir_ind/${assembly}/${batch}/${run}/MT/QC/\${sample_name}_per_base_coverage.tsv ]; then
		per_base_coverage=\$(find $params.outdir_ind/${assembly}/*/${run}/MT/QC/ -name \${sample_name}_per_base_coverage.tsv)
		ln -s \$per_base_coverage .
		touch \${sample_name}_MT_Step1_input_tsv.tsv
	else
		Rscript ../../../modules/shift_back.R $MT_shifted_CollectMetrics \${sample_name}_sorted_chrM_Homo_sapiens_assembly38_collect_wgs_metrics_non_control_region.chrM.interval_list.tsv $path_R_libraries
		mv per_base_coverage.tsv \${sample_name}_per_base_coverage.tsv	

		echo "\${sample_name}\t$params.outdir_ind/${assembly}/${batch}/${run}/MT/QC/\${sample_name}_per_base_coverage.tsv\t\${sample_name}" > \${sample_name}_MT_Step1_input_tsv.tsv
	fi
	"""
}
