// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// MT. Collect coverage and performance metrics for BAM file using CollectWgsMetrics
// This tool collects metrics about the fractions of reads that pass base- and mapping-quality filters as well as coverage (read-depth) levels for WGS analyses.


process Picard_CollectWgsMetrics_MT {
        tag "${bam_MT}"

        publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/Mutect2/QC/", mode: 'copy'

        input :
	file ref_genome_MT_file
        file bam_MT
	file bai_MT
	val assembly
	val batch
	val run

        output :
        file '*_collect_wgs_metrics_MT.txt'

        script :
        """
	gatk CollectWgsMetrics \
        -I ${bam_MT} \
        -O ${bam_MT.simpleName}_collect_wgs_metrics_MT.txt \
        -R ${ref_genome_MT_file}
        """
}
