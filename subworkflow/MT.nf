// Nextflow sub-workflow
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the sub-workflow goal and characteristics :
// Call the MT variants
// Complex pipeline largely inspired from The GATK best practices workflow for Mitochondrial short variant discovery (SNVs + Indels) : https://gatk.broadinstitute.org/hc/en-us/articles/4403870837275-Mitochondrial-short-variant-discovery-SNVs-Indels- and Laricchia et al, 2021, bioRXiv : https://www.biorxiv.org/content/10.1101/2021.07.23.453510v1.full.pdf



// Load the modules for the MT workflow

include { bwa_index; bwa_index as bwa_index_shifted } from "./../modules/bwa_index"
include { Extract_MT_Read } from "./../modules/MT_Extract_MT_Read"
include { MT_SamtoFastq } from "./../modules/MT_SamtoFastq"
include { align_to_MT; align_to_MT as align_to_MT_shifted } from "./../modules/MT_align_to_MT"
include { MarkDuplicates; MarkDuplicates as MarkDuplicates_shifted } from "./../modules/MT_MarkDuplicates"
include { Picard_CollectWgsMetrics_MT; Picard_CollectWgsMetrics_MT as Picard_CollectWgsMetrics_MT_shifted } from "./../modules/MT_Picard_CollectWgsMetrics_MT"
include { MT_call_variants; MT_call_variants as MT_call_variants_shifted } from "./../modules/MT_call_variants"
include { MT_Liftover } from "./../modules/MT_Liftover"
include { MT_MergeVcfs } from "./../modules/MT_MergeVcfs"
include { MT_Merge_stat_file } from "./../modules/MT_Merge_stat_file"
include { MT_Filter_Mutect_Calls } from "./../modules/MT_Filter_Mutect_Calls"
include { MT_LeftAlignAndTrimVariants } from "./../modules/MT_LeftAlignAndTrimVariants"
include { MT_FilterOut_sites } from "./../modules/MT_FilterOut_sites"
include { list_vcfs_txt } from "./../modules/list_vcfs_txt"
include { merge_samples } from "./../modules/merge_samples"

// MT Workflow

workflow MT {

	// Load the parameters and files
	run                                     = params.run
	batch                                   = params.batch
	assembly                                = params.assembly
        assembly_MT                             = params.assembly_MT
	chrM					= params.chrM
	reference                               = file (params.ref)
	reference_index                      	= file (params.ref_index)

	// Load the MT specific files
	Mitochondrial_chromosome                = params.Mitochondrial_chromosome
	ref_MT_fasta                            = file (params.ref_genome_MT)
	ref_MT_fasta_index                      = file (params.ref_genome_MT_index)
	ref_MT_fasta_dict                       = file (params.ref_genome_MT_dict)
	ref_MT_shifted_fasta                    = file (params.ref_genome_MT_shifted)
	ref_MT_shifted_fasta_index              = file (params.ref_genome_MT_shifted_index)
	ref_MT_shifted_fasta_dict               = file (params.ref_genome_MT_shifted_dict)
	blacklist_sites_hg38_MT_file            = file (params.blacklist_sites_hg38_MT)
	blacklist_sites_hg38_MT_index_file      = file (params.blacklist_sites_hg38_MT_index)
	ShiftBack_chain                         = file (params.ShiftBack_chain)
	gnomad_MT_frequ				= file (params.gnomad_MT_frequ)
	MT					= params.MT

        take : 
		bam
		bai
	main:
		bwa_index(ref_MT_fasta)
		bwa_index_shifted(ref_MT_shifted_fasta)
		Extract_MT_Read(bam, bai, Mitochondrial_chromosome)
		MT_SamtoFastq(Extract_MT_Read.out)
		align_to_MT(ref_MT_fasta, bwa_index.out, MT_SamtoFastq.out.fastq_MT)
		align_to_MT_shifted(ref_MT_shifted_fasta, bwa_index_shifted.out, MT_SamtoFastq.out.fastq_MT)
		MarkDuplicates(align_to_MT.out.align_to_MT_bam, align_to_MT.out.align_to_MT_bai)
		MarkDuplicates_shifted(align_to_MT_shifted.out.align_to_MT_bam, align_to_MT_shifted.out.align_to_MT_bai)
		Picard_CollectWgsMetrics_MT(ref_MT_fasta, align_to_MT.out.align_to_MT_bam, align_to_MT.out.align_to_MT_bai, assembly, batch, run)
		Picard_CollectWgsMetrics_MT_shifted(ref_MT_shifted_fasta, align_to_MT_shifted.out.align_to_MT_bam, align_to_MT_shifted.out.align_to_MT_bai, assembly, batch, run)
		MT_call_variants(ref_MT_fasta, ref_MT_fasta_index, ref_MT_fasta_dict, MarkDuplicates.out.bam, MarkDuplicates.out.bai, Mitochondrial_chromosome)
		MT_call_variants_shifted(ref_MT_shifted_fasta, ref_MT_shifted_fasta_index, ref_MT_shifted_fasta_dict, MarkDuplicates_shifted.out.bam, MarkDuplicates_shifted.out.bai, Mitochondrial_chromosome)
		MT_Liftover(MT_call_variants_shifted.out.Mutect2_vcf, MT_call_variants_shifted.out.Mutect2_vcf_index, ref_MT_fasta, ref_MT_fasta_dict, bwa_index.out, ShiftBack_chain, assembly, batch, run)
		MT_MergeVcfs(MT_Liftover.out.lifted_vcf.collect(), MT_call_variants.out.Mutect2_vcf, assembly, batch, run)
		MT_Merge_stat_file(MT_call_variants.out.Mutect2_stat, MT_call_variants_shifted.out.Mutect2_stat.collect())
		MT_Filter_Mutect_Calls(ref_MT_fasta, ref_MT_fasta_index, ref_MT_fasta_dict, MT_MergeVcfs.out.vcf, MT_MergeVcfs.out.index, MT_Merge_stat_file.out.collect())
		MT_LeftAlignAndTrimVariants(ref_MT_fasta, ref_MT_fasta_index, ref_MT_fasta_dict, MT_Filter_Mutect_Calls.out.vcf, MT_Filter_Mutect_Calls.out.index)
		MT_FilterOut_sites(ref_MT_fasta, ref_MT_fasta_index, ref_MT_fasta_dict, MT_LeftAlignAndTrimVariants.out.vcf, MT_LeftAlignAndTrimVariants.out.index, blacklist_sites_hg38_MT_file, blacklist_sites_hg38_MT_index_file, assembly, batch, run)
		list_vcfs_txt(MT_FilterOut_sites.out.collect(), assembly, batch, run, MT)
		merge_samples(list_vcfs_txt.out, assembly, batch, run, MT)
}



// Modules that may be reintegrated if not included in hail 

//include { Bcftools_stats } from "./../modules/Bcftools_stats"
//include { Vcftools_TsTv_by_qual } from "./../modules/Vcftools_TsTv_by_qual"
//include { annotation_table_merged } from "./../modules/annotation_table_merged"
//include { MT_data_organization } from "./../modules/MT_data_organization"
//include { multiqc_pop } from "./../modules/multiqc_pop"

//        vep_cache_merged                        = params.vep_cache_merged
//        vep_cache_merged_version                = params.vep_cache_merged_version
//        CADD_1_6_whole_genome_SNVs              = file (params.CADD_1_6_whole_genome_SNVs)
//        CADD_1_6_whole_genome_SNVs_index        = file (params.CADD_1_6_whole_genome_SNVs_index)
//        CADD_1_6_InDels                         = file (params.CADD_1_6_InDels)
//        CADD_1_6_InDels_index                   = file (params.CADD_1_6_InDels_index)

//              QC1 = Bcftools_stats(merge_samples.out.vcf, merge_samples.out.index, assembly, run)
//              QC2 = Vcftools_TsTv_by_qual(merge_samples.out.vcf, merge_samples.out.index, assembly, run)
//              annotation_table_merged(merge_samples.out.vcf, merge_samples.out.index, vep_cache_merged, vep_cache_merged_version, assembly, run, assembly_MT, CADD_1_6_whole_genome_SNVs, CADD_1_6_whole_genome_SNVs_index, CADD_1_6_InDels, CADD_1_6_InDels_index, chrM, MT)
//              MT_data_organization(gnomad_MT_frequ, merge_samples.out.vcf, annotation_table_merged.out.annot_table_merged_R, assembly, run)
//              quality_metrics = QC1.concat(QC2, annotation_table_merged.out.vep_merged_stat).collect()
//              quality_metrics = QC1.concat(QC2).collect()
//              multiqc_pop(quality_metrics, assembly, run, Mitochondrial_chromosome)
