// Nextflow sub-workflow
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the sub-workflow goal and characteristics :
// Index the reference genome (specified in the nextflow.config file)
// Align, sort and index the fastq for each sample (fastq --> bam)


// Load the modules for the ALN workflow
include { align_sort_output_bam } from "./../modules/align_sort_output_bam"
include { bwa_index; bwa_index as bwa_index_shifted } from "./../modules/bwa_index"

// ALN workflow
workflow ALN {

	// Load the parameters and files
	run             = params.run
	batch           = params.batch
	assembly        = params.assembly
	reference       = file (params.ref)

	Channel
    		.fromFilePairs(params.reads )
    		.set {read_pairs_ch}

	main:
		bwa_index(reference)
		align_sort_output_bam(reference, bwa_index.out, read_pairs_ch, assembly, batch, run)
	emit :
		reference_index 	= bwa_index.out.collect()
		bam_sorted 		= align_sort_output_bam.out.samples_bam
		bam_sorted_index 	= align_sort_output_bam.out.samples_bam_index
}

