// Nextflow sub-workflow
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the sub-workflow goal and characteristics :
// Call the SNV variants


// Load the modules for the SNV workflow

//include { gnomad_frequency_table } from "./../modules/gnomad_frequency_table"
include { deepvariant_call } from "./../modules/deepvariant"
include { gvcfs_txt } from "./../modules/gvcf_txt"
include { GLnexus_cli } from "./../modules/GLnexus"
include { bcf_to_vcf } from "./../modules/bcf_to_vcf"
//include { split_vcf_by_chr } from "./../modules/split_vcf_by_chr"
include { Bcftools_stats } from "./../modules/Bcftools_stats"
include { Vcftools_TsTv_by_qual } from "./../modules/Vcftools_TsTv_by_qual"
//include { annotation_table_merged } from "./../modules/annotation_table_merged"
//include { SNV_data_organization } from "./../modules/SNV_data_organization"
include { multiqc_pop } from "./../modules/multiqc_pop"
include { plink_sex_inference } from "./../modules/plink_sex_inference"
include { sample_QC } from "./../modules/sample_QC"
include { count_variants_vcftools } from "./../modules/count_variants_vcftools"
include { count_variants_gatk } from "./../modules/count_variants_gatk"
include { count_variants_gatk_2 } from "./../modules/count_variants_gatk_2"
include { count_bcftools_stats } from "./../modules/count_bcftools_stats"

include { Hail_sample_QC } from "./../modules/Hail_sample_QC"
include { Hail_variant_QC } from "./../modules/Hail_variant_QC"


// SNV workflow

workflow SNV {

	// Load the parameters and files
	run             			= params.run
	batch           			= params.batch
	assembly        			= params.assembly
	assembly_hg				= params.assembly_hg
	reference       			= file (params.ref)
	reference_index 			= file (params.ref_index)
	vep_cache_merged 			= file (params.vep_cache_merged)
	vep_cache_merged_version 		= params.vep_cache_merged_version
	CADD_1_6_whole_genome_SNVs 		= file (params.CADD_1_6_whole_genome_SNVs)
	CADD_1_6_whole_genome_SNVs_index 	= file (params.CADD_1_6_whole_genome_SNVs_index)
	CADD_1_6_InDels 			= file (params.CADD_1_6_InDels)
	CADD_1_6_InDels_index 			= file (params.CADD_1_6_InDels_index)
	SNV					= params.SNV
	gnomad_SNV_vcf				= file (params.gnomad_SNV_vcf)
	gnomad_SNV_index			= file (params.gnomad_SNV_index)
	chr					= params.chrom

	// Workflow start here
	take : 
		bam
		bai
		mosdepth

	main :
//		gnomad_frequency_table(gnomad_SNV_vcf, gnomad_SNV_index, chr)

		deepvariant_call(reference, reference_index, bam, bai, assembly, batch, run)
		gvcfs_txt(deepvariant_call.out.deepvariant_gvcf.collect(), assembly, batch, run)
		GLnexus_cli(gvcfs_txt.out, run)
		bcf_to_vcf(GLnexus_cli.out, assembly, batch, run)

		count_variants_gatk(deepvariant_call.out.deepvariant_vcf, deepvariant_call.out.deepvariant_vcf_index, assembly, batch, run)
                count_variants_gatk_2(count_variants_gatk.out.collect(), assembly, batch, run)
                count_variants_vcftools(bcf_to_vcf.out.vcf, bcf_to_vcf.out.index, assembly, batch, run)
                count_bcftools_stats(bcf_to_vcf.out.vcf, bcf_to_vcf.out.index, assembly, batch, run)
                plink_sex_inference(bcf_to_vcf.out.vcf, assembly_hg, assembly, batch, run)
		sample_QC(plink_sex_inference.out, assembly, batch, run, count_variants_vcftools.out, count_variants_gatk_2.out, count_bcftools_stats.out, mosdepth)

//		split_vcf_by_chr(bcf_to_vcf.out.vcf, assembly, batch, run, chr)
//                annotation_table_merged(bcf_to_vcf.out.vcf, bcf_to_vcf.out.index, vep_cache_merged, vep_cache_merged_version, assembly, run, assembly, CADD_1_6_whole_genome_SNVs, CADD_1_6_whole_genome_SNVs_index, CADD_1_6_InDels, CADD_1_6_InDels_index, chr, SNV)
//                SNV_data_organization(gnomad_frequency_table.out.collect(), split_vcf_by_chr.out.vcf_onechr, annotation_table_merged.out.annot_table_merged_R.collect(), assembly, run, sample_QC.out.sample_QC_file)

                QC1 = Bcftools_stats(bcf_to_vcf.out.vcf, bcf_to_vcf.out.index, assembly, run)
                QC2 = Vcftools_TsTv_by_qual(bcf_to_vcf.out.vcf, bcf_to_vcf.out.index, assembly, run)
//                quality_metrics = QC1.concat(QC2, annotation_table_merged.out.vep_merged_stat).collect()
                quality_metrics = QC1.concat(QC2).collect()
                multiqc_pop(quality_metrics, assembly, run, SNV)

		Hail_sample_QC(bcf_to_vcf.out.vcf, assembly, batch, run)
                Hail_variant_QC(Hail_sample_QC.out.vcf_sample_filtered, assembly, batch, run)


	emit :
		sample_sex_file=sample_QC.out.sample_QC_file
}

