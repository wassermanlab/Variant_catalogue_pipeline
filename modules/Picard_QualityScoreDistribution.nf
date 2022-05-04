// Quality check of the bam file with Picard QualityScoreDistribution
// This tool is used for determining the overall 'quality' for a library in a given run. To that effect, it outputs a chart and tables indicating the range of quality scores and the total numbers of bases corresponding to those scores. 
// (!) R is necessary for the chart


process Picard_QualityScoreDistribution {
        tag "${bam.simpleName}"
 
	publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/QC/Individuals/${bam.simpleName}/Picard_Metrics/", mode: 'copy'
	
	input :
	file bam
	file bai
	val assembly
	val batch
	val run

	output :
	file '*_qual_score_dist.*' 

	script :
	"""
	picard "-Xmx2G" QualityScoreDistribution \
	I=${bam} \
	O=${bam.simpleName}_qual_score_dist.txt \
	CHART= ${bam.simpleName}_qual_score_dist.pdf
	"""
}

