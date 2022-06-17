// Nextflow sub-workflow
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the sub-workflow goal and characteristics :
// Hail producing several graphs and filtering outliers samples and variants

// Load the modules for the SNV workflow

include { Hail_sample_QC } from "./../modules/Hail_sample_QC"
include { Hail_variant_QC } from "./../modules/Hail_variant_QC"

include { annotation_table_merged as SNV_annotation_table_merged; annotation_table_merged as MT_annotation_table_merged} from "./../modules/annotation_table_merged"

include { split_tsv_by_chr } from "./../modules/split_tsv_by_chr"
include { gnomad_frequency_table } from "./../modules/gnomad_frequency_table"
include { SNV_data_organization } from "./../modules/SNV_data_organization"

//include { Hail_variant_MT_QC } from "./../modules/Hail_variant_MT_QC"
//include { MT_data_organization } from "./../modules/MT_data_organization"

// Hail workflow

workflow Hail {

	// Load the parameters and files
	run             			= params.run
	batch           			= params.batch
	assembly        			= params.assembly
	assembly_hg				= params.assembly_hg
        chr                                     = params.chrom
        vep_cache_merged                        = file (params.vep_cache_merged)
        vep_cache_merged_version                = params.vep_cache_merged_version
        CADD_1_6_whole_genome_SNVs              = file (params.CADD_1_6_whole_genome_SNVs)
        CADD_1_6_whole_genome_SNVs_index        = file (params.CADD_1_6_whole_genome_SNVs_index)
        CADD_1_6_InDels                         = file (params.CADD_1_6_InDels)
        CADD_1_6_InDels_index                   = file (params.CADD_1_6_InDels_index)
	spliceai_snv				= file (params.spliceai_snv)
	spliceai_snv_index			= file (params.spliceai_snv_index)
	spliceai_indel				= file (params.spliceai_indel)
	spliceai_indel_index			= file (params.spliceai_indel_index)
        SNV                                     = params.SNV
	gnomad_SNV_vcf                          = file (params.gnomad_SNV_vcf)
	gnomad_SNV_index                        = file (params.gnomad_SNV_index)
        severity_table                          = file (params.severity_table)

	// Workflow start here
	take : 
		SNV_vcf
		sample_sex_file
	main :
		Hail_sample_QC(SNV_vcf, assembly, batch, run)
                Hail_variant_QC(Hail_sample_QC.out.vcf_sample_filtered, sample_sex_file, assembly, batch, run)

		gnomad_frequency_table(gnomad_SNV_vcf, gnomad_SNV_index, chr)
		split_tsv_by_chr(Hail_variant_QC.out.SNV_frequ_tot_xx_xy_tsv, assembly, batch, run)
		
                SNV_annotation_table_merged(Hail_variant_QC.out.SNV_vcf_samples_variants_filtered, Hail_variant_QC.out.SNV_index, vep_cache_merged, vep_cache_merged_version, assembly, run, assembly, CADD_1_6_whole_genome_SNVs, CADD_1_6_whole_genome_SNVs_index, CADD_1_6_InDels, CADD_1_6_InDels_index, spliceai_snv, spliceai_snv_index, spliceai_indel, spliceai_indel_index, chr, SNV)
		SNV_data_organization(gnomad_frequency_table.out.collect(), split_tsv_by_chr.out.collect(), SNV_annotation_table_merged.out.annot_table_merged_R.collect(), assembly, run, chr, Hail_variant_QC.out.SNV_vcf_samples_variants_filtered, severity_table)
}

// Solenne : To have independant processes and allow users of the pipeline to use the MT, the SV and the SNV pipelines separately, the MT data filtering steps using hail are merged with the MT pipeline
// To discuss in the future : Should the SNV and SV data filtering be separate and included into each pipeline (SNV < 50bp and SV >= 50 bp to avoid overlaps)?

//include { Hail_variant_MT_QC } from "./../modules/Hail_variant_MT_QC"
//include { MT_data_organization } from "./../modules/MT_data_organization"

//        gnomad_MT_frequ                         = file (params.gnomad_MT_frequ)
//        MT                                      = params.MT
//        chrM                                    = params.chrM


//                Hail_variant_MT_QC(Hail_sample_QC.out.vcf_sample_filtered, MT_vcf, sample_sex_file, assembly, batch, run)
//                MT_annotation_table_merged(Hail_variant_MT_QC.out.MT_vcf_samples_variants_filtered, Hail_variant_MT_QC.out.MT_index, vep_cache_merged, vep_cache_merged_version, assembly, run, assembly, CADD_1_6_whole_genome_SNVs, CADD_1_6_whole_genome_SNVs_index, CADD_1_6_InDels, CADD_1_6_InDels_index, chrM, MT)
//                MT_data_organization(gnomad_MT_frequ, Hail_variant_MT_QC.out.MT_frequ_tot_xx_xy_tsv, annotation_table_merged.out.annot_table_merged_R, assembly, run)                          
