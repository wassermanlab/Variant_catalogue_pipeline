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
include { MarkDuplicates_index; MarkDuplicates_index as MarkDuplicates_index_shifted } from "./../modules/MT_MarkDuplicates_index"

include { Picard_CollectWgsMetrics_MT; Picard_CollectWgsMetrics_MT as Picard_CollectWgsMetrics_MT_shifted } from "./../modules/MT_Picard_CollectWgsMetrics_MT"
include { shift_back } from "./../modules/shift_back"
include { MT_call_variants; MT_call_variants as MT_call_variants_shifted } from "./../modules/MT_call_variants"
include { MT_Liftover } from "./../modules/MT_Liftover"
include { MT_MergeVcfs } from "./../modules/MT_MergeVcfs"
include { MT_norm } from "./../modules/MT_norm"
include { MT_Merge_stat_file } from "./../modules/MT_Merge_stat_file"
include { MT_Filter_Mutect_Calls } from "./../modules/MT_Filter_Mutect_Calls"
include { MT_LeftAlignAndTrimVariants } from "./../modules/MT_LeftAlignAndTrimVariants"
include { MT_FilterOut_sites } from "./../modules/MT_FilterOut_sites"
include { MT_haplocheck } from "./../modules/MT_haplocheck"
include { MT_Step1_input_tsv } from "./../modules/MT_Step1_input_tsv"
include { MT_Step2_participant_data } from "./../modules/MT_Step2_participant_data"
include { MT_Step3_metadata } from "./../modules/MT_Step3_metadata"
include { MT_Step3_metadata_sample } from "./../modules/MT_Step3_metadata_sample"
include { Hail_variant_MT_QC } from "./../modules/Hail_variant_MT_QC"
include { annotation_table_merged } from "./../modules/annotation_table_merged"

include { MT_data_organization } from "./../modules/MT_data_organization"


// MT Workflow

workflow MT {

	// Load the parameters and files
	run                                     = params.run
	batch                                   = params.batch
	assembly                                = params.assembly
        assembly_MT                             = params.assembly_MT
	chrM					= params.chrM
	reference                               = file (params.ref)
	vep_cache_merged                        = params.vep_cache_merged
	dir_plugin				= params.dir_plugin
        vep_cache_merged_version                = params.vep_cache_merged_version
        CADD_1_6_whole_genome_SNVs              = file (params.CADD_1_6_whole_genome_SNVs)
        CADD_1_6_whole_genome_SNVs_index        = file (params.CADD_1_6_whole_genome_SNVs_index)
        CADD_1_6_InDels                         = file (params.CADD_1_6_InDels)
        CADD_1_6_InDels_index                   = file (params.CADD_1_6_InDels_index)
	spliceai_snv				= file (params.spliceai_snv)
	spliceai_snv_index			= file (params.spliceai_snv_index)
	spliceai_indel				= file (params.spliceai_indel)
	spliceai_indel_index			= file (params.spliceai_indel_index)
	severity_table				= file (params.severity_table)
	path_R_libraries			= params.path_R_libraries
	haplocheck_path				= params.haplocheck_path

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
	non_control_region_interval_list	= file (params.non_control_region_interval_list)
	control_region_shifted_reference_interval_list = file (params.control_region_shifted_reference_interval_list)
	gnomad_MT_frequ				= file (params.gnomad_MT_frequ)
	MT					= params.MT
	pon_predictions_table			= file (params.pon_predictions_table)
	artifact_prone_sites_bed		= file (params.artifact_prone_sites_bed)
        mitotip_predictions_table               = file (params.mitotip_predictions_table)

//A copy of the reference was done locally because Hail could not read the copy located in a common directory
	GRCh38_MT_local_fasta			= file (params.GRCh38_MT_local_fasta)
	GRCh38_MT_local_fai			= file (params.GRCh38_MT_local_fai)
	

        take : 
		bam
		bai
		mosdepth
	main:
		bwa_index(ref_MT_fasta)
		bwa_index_shifted(ref_MT_shifted_fasta)

                // Sample specific (Do not need to be run for a previously processed sample)
		Extract_MT_Read(bam, bai, Mitochondrial_chromosome, assembly, batch, run)
		MT_SamtoFastq(Extract_MT_Read.out, assembly, batch, run)
		//
                align_to_MT(ref_MT_fasta, bwa_index.out, MT_SamtoFastq.out.fastq_MT, 
                assembly, batch, run)
		//
                align_to_MT_shifted(ref_MT_shifted_fasta, bwa_index_shifted.out,
                MT_SamtoFastq.out.fastq_MT, assembly, batch, run)
		//
                MarkDuplicates(align_to_MT.out.align_to_MT_bam, 
                align_to_MT.out.align_to_MT_bai, assembly, batch, run)
                //
                MarkDuplicates_index(align_to_MT.out.align_to_MT_bam, 
                align_to_MT.out.align_to_MT_bai, assembly, batch, run, MarkDuplicates.out.bam.collect())
                //    
		MarkDuplicates_shifted(align_to_MT_shifted.out.align_to_MT_bam,
                align_to_MT_shifted.out.align_to_MT_bai, assembly, batch, run)
		//
                MarkDuplicates_index_shifted(align_to_MT_shifted.out.align_to_MT_bam,
                align_to_MT_shifted.out.align_to_MT_bai, assembly, batch, run,
                MarkDuplicates_shifted.out.bam.collect())
		//
                Picard_CollectWgsMetrics_MT(ref_MT_fasta, ref_MT_fasta_index,
                non_control_region_interval_list, align_to_MT.out.align_to_MT_bam,
                align_to_MT.out.align_to_MT_bai, assembly, batch, run)
		//
                Picard_CollectWgsMetrics_MT_shifted(ref_MT_shifted_fasta,
                ref_MT_shifted_fasta_index, control_region_shifted_reference_interval_list,
                align_to_MT_shifted.out.align_to_MT_bam, align_to_MT_shifted.out.align_to_MT_bai,
                assembly, batch, run)
		//
                shift_back(Picard_CollectWgsMetrics_MT_shifted.out,
                Picard_CollectWgsMetrics_MT.out.collect(), assembly, batch, run, path_R_libraries)
		MT_Step1_input_tsv(shift_back.out.Sample_MT_Step1_input_tsv.collect(), assembly, batch, run)
		MT_call_variants(ref_MT_fasta, ref_MT_fasta_index, ref_MT_fasta_dict,
                MarkDuplicates.out.bam, MarkDuplicates_index.out.bai.collect(),
                Mitochondrial_chromosome, assembly, batch, run)
		MT_call_variants_shifted(ref_MT_shifted_fasta, ref_MT_shifted_fasta_index, 
                ref_MT_shifted_fasta_dict, MarkDuplicates_shifted.out.bam,
                MarkDuplicates_index_shifted.out.bai.collect(), Mitochondrial_chromosome, assembly, batch, run)
		MT_Liftover(MT_call_variants_shifted.out.Mutect2_vcf,
                MT_call_variants_shifted.out.Mutect2_vcf_index, ref_MT_fasta, ref_MT_fasta_dict,
                bwa_index.out, ShiftBack_chain, assembly, batch, run)
		MT_MergeVcfs(MT_Liftover.out.lifted_vcf.collect(), MT_call_variants.out.Mutect2_vcf, assembly, batch, run)
		MT_norm(MT_MergeVcfs.out.vcf, assembly, batch, run)
		//
                MT_Merge_stat_file(MT_call_variants.out.Mutect2_stat,
                //
                MT_call_variants_shifted.out.Mutect2_stat.collect(), assembly, batch, run)
		//
                MT_Filter_Mutect_Calls(ref_MT_fasta, ref_MT_fasta_index, ref_MT_fasta_dict,
                MT_norm.out.vcf, MT_norm.out.index, MT_Merge_stat_file.out.collect(), assembly, batch, run)
		//
                MT_LeftAlignAndTrimVariants(ref_MT_fasta, ref_MT_fasta_index, ref_MT_fasta_dict, 
                //
                MT_Filter_Mutect_Calls.out.vcf, MT_Filter_Mutect_Calls.out.index, assembly, batch, run)
		//
                MT_FilterOut_sites(ref_MT_fasta, ref_MT_fasta_index, ref_MT_fasta_dict,
                MT_LeftAlignAndTrimVariants.out.vcf, MT_LeftAlignAndTrimVariants.out.index,
                blacklist_sites_hg38_MT_file, blacklist_sites_hg38_MT_index_file, assembly, batch, run)
                //
                MT_haplocheck(MT_FilterOut_sites.out.vcf, assembly, batch, run, haplocheck_path)
                MT_Step3_metadata_sample(mosdepth, MT_haplocheck.out.file, assembly, batch,
                run, path_R_libraries)

                // Aggregated steps (Need to be run everytime a new sample is added to the cohort)
		MT_Step2_participant_data(MT_FilterOut_sites.out.sample_MT_Step2_participant_data.collect(), MT_FilterOut_sites.out.Sample_list.collect(), assembly, batch, run)
		MT_Step3_metadata(MT_Step3_metadata_sample.out.collect(), assembly, batch, run)
		Hail_variant_MT_QC(MT_Step1_input_tsv.out, MT_Step2_participant_data.out.MT_Step2_participant_data_tsv, MT_Step2_participant_data.out.participants_to_subset_txt, MT_Step3_metadata.out, assembly, batch, run, pon_predictions_table, artifact_prone_sites_bed, GRCh38_MT_local_fasta, GRCh38_MT_local_fai, mitotip_predictions_table)
                annotation_table_merged(Hail_variant_MT_QC.out.Hail_MT_frequ_only, Hail_variant_MT_QC.out.vcf_index, vep_cache_merged, vep_cache_merged_version, assembly, run, assembly_MT, CADD_1_6_whole_genome_SNVs, CADD_1_6_whole_genome_SNVs_index, CADD_1_6_InDels, CADD_1_6_InDels_index, spliceai_snv, spliceai_snv_index, spliceai_indel, spliceai_indel_index, chrM, MT, reference, dir_plugin)

                MT_data_organization(gnomad_MT_frequ, annotation_table_merged.out.annotation_vcf, assembly, run, severity_table)
}
