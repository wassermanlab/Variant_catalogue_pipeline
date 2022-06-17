// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// Quality check of the bam file with Picard CollectWgsMetrics
// This tool collects metrics about the fractions of reads that pass base- and mapping-quality filters as well as coverage (read-depth) levels for WGS analyses.


process Picard_CollectWgsMetrics {
        tag "${bam.simpleName}"
 
	publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/QC/Individuals/${bam.simpleName}/Picard_Metrics/", mode: 'copyNoFollow'

	input :
	file bam
	file bai
	path ref_genome 
	path reference_index
	val assembly
	val batch
	val run

	output :
	file '*_collect_wgs_metrics.txt' 

	script :
	"""
	if [ -a $params.outdir_ind/${assembly}/*/${run}/QC/Individuals/${bam.simpleName}/Picard_Metrics/${bam.simpleName}_collect_wgs_metrics.txt ]; then
		picard_collect_wgs_metrics=\$(find $params.outdir_ind/${assembly}/*/${run}/QC/Individuals/${bam.simpleName}/Picard_Metrics/ -name ${bam.simpleName}_collect_wgs_metrics.txt)
		ln -s \$picard_collect_wgs_metrics .
	else
		gatk CollectWgsMetrics \
		--java-options "-Xmx8G" \
		-I ${bam} \
		-O ${bam.simpleName}_collect_wgs_metrics.txt \
		-R ${ref_genome}
	fi
	"""
}
