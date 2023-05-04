// Nextflow sub-workflow
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the sub-workflow goal and characteristics :
// Index the reference genome (specified in the nextflow.config file)
// Align, sort and index the fastq for each sample (fastq --> bam)
// Genereate Quality control (QC) metrics for each sample using different software
// The results are aggregated using multiQC


// Load the modules for the mapping workflow
include { align_sort_output_bam } from "./../modules/align_sort_output_bam"
include { bwa_index } from "./../modules/bwa_index"

include {fastqc} from "./../modules/fastqc"
include {Mosdepth} from "./../modules/mosdepth"
include {Picard_CollectWgsMetrics} from "./../modules/Picard_CollectWgsMetrics"
include {Picard_CollectAlignmentSummaryMetrics} from "./../modules/Picard_CollectAlignmentSummaryMetrics"
include {Picard_QualityScoreDistribution} from "./../modules/Picard_QualityScoreDistribution"
include {multiqc_indiv} from "./../modules/multiqc_indiv"


// mapping workflow
workflow Mapping {

	// Load the parameters and files
        batch           = params.batch
        assembly        = params.assembly
        run             = params.run
        outdir_ind      = file (params.outdir_ind)
        reference       = file (params.ref)
        reference_index = file (params.ref_index)

	Channel
    		.fromFilePairs(params.reads )
    		.set {read_pairs_ch}

	main:
		bwa_index(reference)
		align_sort_output_bam(reference, bwa_index.out, read_pairs_ch, assembly, batch, run)

                q1              = fastqc(read_pairs_ch, outdir_ind, assembly, batch, run)
                q2              = Mosdepth(align_sort_output_bam.out.samples_bam, align_sort_output_bam.out.samples_bam_index, assembly, batch, run)
                q3              = Picard_CollectWgsMetrics(align_sort_output_bam.out.samples_bam, align_sort_output_bam.out.samples_bam_index, reference, reference_index, assembly, batch, run)
                q4              = Picard_CollectAlignmentSummaryMetrics(align_sort_output_bam.out.samples_bam, align_sort_output_bam.out.samples_bam_index, assembly, batch, run)
                q5              = Picard_QualityScoreDistribution(align_sort_output_bam.out.samples_bam, align_sort_output_bam.out.samples_bam_index, assembly, batch, run)
                quality_metrics = q1.concat(q2.all_files,q3,q4,q5).collect()
                multiqc_indiv(quality_metrics, assembly, batch, run)

	emit :
		reference_index 	= bwa_index.out.collect()
		bam_sorted 		= align_sort_output_bam.out.samples_bam
		bam_sorted_index 	= align_sort_output_bam.out.samples_bam_index
                mosdepth_output 	= Mosdepth.out.summary_stat.collect()

}


