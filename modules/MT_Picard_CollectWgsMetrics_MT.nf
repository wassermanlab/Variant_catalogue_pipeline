// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// MT. Collect coverage and performance metrics for BAM file using CollectWgsMetrics
// This tool collects metrics about the fractions of reads that pass base- and mapping-quality filters as well as coverage (read-depth) levels for WGS analyses.


process Picard_CollectWgsMetrics_MT {
        tag "${bam_MT}"

        publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/MT/QC/${bam_MT.simpleName}/", mode: 'copyNoFollow'

        input :
	file ref_genome_MT_file
        file ref_genome_MT_file_index
        file interval_list
	file bam_MT
	file bai_MT
	val assembly
	val batch
	val run

        output :
        file '*.tsv'

        script :
        """
	sample_name=\$(echo ${bam_MT.simpleName} | cut -d _ -f 1)
	if [ -a $params.outdir_ind/${assembly}/*/${run}/MT/QC/*/\${sample_name}*_collect_wgs_metrics_${interval_list}.tsv ]; then
		metric=\$(find $params.outdir_ind/${assembly}/*/${run}/MT/QC/*/ -name \${sample_name}*_collect_wgs_metrics_${interval_list}.tsv)
		ln -s \$metric .
	else
		gatk CollectHsMetrics \
        	--java-options "-Xmx8G" \
		-I ${bam_MT} \
        	--PER_BASE_COVERAGE ${bam_MT.simpleName}_collect_wgs_metrics_${interval_list}.tsv \
        	-R ${ref_genome_MT_file} \
		-O ${bam_MT.simpleName}.metrics \
		-TI $interval_list \
		-BI $interval_list \
		--SAMPLE_SIZE 1 
	fi
	"""
}
