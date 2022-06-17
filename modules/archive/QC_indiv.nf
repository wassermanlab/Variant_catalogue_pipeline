// Nextflow sub-workflow
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the sub-workflow goal and characteristics :
// Genereate Quality control (QC) metrics for each sample using different software
// The results are aggregated using multiQC


// Load the modules for the MT workflow

include {fastqc} from "./../modules/fastqc"
include {Mosdepth} from "./../modules/mosdepth"
include {Picard_CollectWgsMetrics} from "./../modules/Picard_CollectWgsMetrics"
include {Picard_CollectAlignmentSummaryMetrics} from "./../modules/Picard_CollectAlignmentSummaryMetrics"
include {Picard_QualityScoreDistribution} from "./../modules/Picard_QualityScoreDistribution"
include {multiqc_indiv} from "./../modules/multiqc_indiv"


workflow QC_indiv {

	// Load the parameters and files
	samples 	= Channel
		.fromFilePairs (params.reads)
		.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }    

	batch 		= params.batch
	assembly        = params.assembly
	run		= params.run
	outdir_ind 	= file (params.outdir_ind)
	reference	= file (params.ref)
	reference_index	= file (params.ref_index)

	// Workflow start here
	take :	
		bam
		bai
	main:
		q1 		= fastqc(samples, outdir_ind, assembly, batch, run)
		q2 		= Mosdepth(bam, bai, assembly, batch, run) 
		q3 		= Picard_CollectWgsMetrics(bam, bai, reference, reference_index, assembly, batch, run)
		q4 		= Picard_CollectAlignmentSummaryMetrics(bam, bai, assembly, batch, run)
		q5 		= Picard_QualityScoreDistribution(bam, bai, assembly, batch, run)
		quality_metrics	= q1.concat(q2.all_files,q3,q4,q5).collect()
		multiqc_indiv(quality_metrics, assembly, batch, run)
	emit:
		mosdepth_output = Mosdepth.out.summary_stat.collect()
}

