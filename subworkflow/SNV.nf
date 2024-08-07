// Nextflow sub-workflow
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developed to build the IBVL, a background variant library

// Overview of the sub-workflow goal and characteristics :
// Call the SNV variants
// Include some quality controls (QC) steps
// Hail producing several graphs and filtering outliers samples and variants

// Load the modules for the SNV workflow

include { deepvariant_call } from "./../modules/deepvariant"
include { list_vcfs_txt } from "./../modules/list_vcfs_txt"
include { GLnexus_cli } from "./../modules/GLnexus"
include { bcf_to_vcf } from "./../modules/bcf_to_vcf"

include { Hail_sample_QC } from "./../modules/Hail_sample_QC"
include { Hail_variant_QC } from "./../modules/Hail_variant_QC"

include { annotation_table_merged as SNV_annotation_table_merged; annotation_table_merged as MT_annotation_table_merged} from "./../modules/annotation_table_merged"

include { SNV_data_organization } from "./../modules/SNV_data_organization"

// SNV workflow

workflow SNV {

	// Load the parameters and files
	run             			= params.run
	batch           			= params.batch
	assembly        			= params.assembly
	reference       			= file (params.ref)
	reference_index 			= file (params.ref_index)
        SNV                                     = params.SNV

        chr                                     = params.chrom
        vep_cache_merged                        = file (params.vep_cache_merged)
        vep_cache_merged_version                = params.vep_cache_merged_version
	dir_plugin				= file (params.dir_plugin)
        CADD_1_6_whole_genome_SNVs              = file (params.CADD_1_6_whole_genome_SNVs)
        CADD_1_6_whole_genome_SNVs_index        = file (params.CADD_1_6_whole_genome_SNVs_index)
        CADD_1_6_InDels                         = file (params.CADD_1_6_InDels)
        CADD_1_6_InDels_index                   = file (params.CADD_1_6_InDels_index)
        spliceai_snv                            = file (params.spliceai_snv)
        spliceai_snv_index                      = file (params.spliceai_snv_index)
        spliceai_indel                          = file (params.spliceai_indel)
        spliceai_indel_index                    = file (params.spliceai_indel_index)
        severity_table                          = file (params.severity_table)
	gnomad_SNV_frequ			= file (params.gnomad_SNV_frequ)

	// Workflow start here
	take : 
		bam
		bai

	main :
		// Sample specific (Do not need to be run for a previously processed sample)
		deepvariant_call(reference, reference_index, bam, bai, assembly, batch, run)

		// Aggregated steps (Need to be run everytime a new sample is added to the cohort)
    //	list_vcfs_txt(deepvariant_call.out.deepvariant_gvcf.collect(), assembly, batch, run, SNV)
	//	GLnexus_cli(list_vcfs_txt.out, run)
	//	bcf_to_vcf(GLnexus_cli.out, assembly, batch, run)

    //            Hail_sample_QC(bcf_to_vcf.out.vcf, assembly,reference,reference_index,batch,run)
    //            Hail_variant_QC(Hail_sample_QC.out.vcf_sample_filtered, Hail_sample_QC.out.filtered_sample_sex, assembly, reference, reference_index,batch,run,chr)
    //            SNV_annotation_table_merged(Hail_variant_QC.out.vcf_SNV_filtered_frequ_only, Hail_variant_QC.out.index_SNV_filtered_frequ_only, vep_cache_merged, vep_cache_merged_version, assembly, run, assembly, CADD_1_6_whole_genome_SNVs, CADD_1_6_whole_genome_SNVs_index, CADD_1_6_InDels, CADD_1_6_InDels_index, spliceai_snv, spliceai_snv_index, spliceai_indel, spliceai_indel_index, SNV, reference, dir_plugin)

	//	SNV_data_organization(gnomad_SNV_frequ, SNV_annotation_table_merged.out.annotation_vcf, assembly, run, severity_table)

	//emit :
	//	sample_sex_file=Hail_sample_QC.out.filtered_sample_sex

}
