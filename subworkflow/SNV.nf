// Nextflow sub-workflow
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the sub-workflow goal and characteristics :
// Call the SNV variants
// Include some quality controls (QC) steps
//	- Plink which defines the sex of each sample based on seevral variables

// Load the modules for the SNV workflow

include { deepvariant_call } from "./../modules/deepvariant"
include { list_vcfs_txt } from "./../modules/list_vcfs_txt"
include { GLnexus_cli } from "./../modules/GLnexus"
include { bcf_to_vcf } from "./../modules/bcf_to_vcf"

include { plink_sex_inference } from "./../modules/plink_sex_inference"
include { sample_QC } from "./../modules/sample_QC"

// SNV workflow

workflow SNV {

	// Load the parameters and files
	run             			= params.run
	batch           			= params.batch
	assembly        			= params.assembly
	assembly_hg				= params.assembly_hg
	reference       			= file (params.ref)
	reference_index 			= file (params.ref_index)
        SNV                                     = params.SNV

	// Workflow start here
	take : 
		bam
		bai
		mosdepth

	main :
		deepvariant_call(reference, reference_index, bam, bai, assembly, batch, run)
		list_vcfs_txt(deepvariant_call.out.deepvariant_gvcf.collect(), assembly, batch, run, SNV)
		GLnexus_cli(list_vcfs_txt.out, run)
		bcf_to_vcf(GLnexus_cli.out, assembly, batch, run)

                plink_sex_inference(bcf_to_vcf.out.vcf, assembly_hg, assembly, batch, run)
		sample_QC(plink_sex_inference.out, assembly, batch, run, mosdepth)

	emit :
		sample_sex_file=sample_QC.out.sample_QC_file
		SNV_vcf = bcf_to_vcf.out.vcf
}


// Steps removed as they are now integrated in Hail
//include { gnomad_frequency_table } from "./../modules/gnomad_frequency_table"
//include { count_variants_vcftools } from "./../modules/count_variants_vcftools"
//include { count_variants_gatk } from "./../modules/count_variants_gatk"
//include { count_variants_gatk_2 } from "./../modules/count_variants_gatk_2"
//include { count_bcftools_stats } from "./../modules/count_bcftools_stats"
//include { split_vcf_by_chr } from "./../modules/split_vcf_by_chr"
//include { Bcftools_stats } from "./../modules/Bcftools_stats"
//include { Vcftools_TsTv_by_qual } from "./../modules/Vcftools_TsTv_by_qual"

//include { SNV_data_organization } from "./../modules/SNV_data_organization"
//include { multiqc_pop } from "./../modules/multiqc_pop"

//        gnomad_SNV_vcf                          = file (params.gnomad_SNV_vcf)
//        gnomad_SNV_index                        = file (params.gnomad_SNV_index)

//              count_variants_gatk(deepvariant_call.out.deepvariant_vcf, deepvariant_call.out.deepvariant_vcf_index, assembly, batch, run)
//                count_variants_gatk_2(count_variants_gatk.out.collect(), assembly, batch, run)
//                count_variants_vcftools(bcf_to_vcf.out.vcf, bcf_to_vcf.out.index, assembly, batch, run)
//                count_bcftools_stats(bcf_to_vcf.out.vcf, bcf_to_vcf.out.index, assembly, batch, run)
//                sample_QC(plink_sex_inference.out, assembly, batch, run, count_variants_vcftools.out, count_variants_gatk_2.out, count_bcftools_stats.out, mosdepth)
//                QC1 = Bcftools_stats(bcf_to_vcf.out.vcf, bcf_to_vcf.out.index, assembly, run)
//                QC2 = Vcftools_TsTv_by_qual(bcf_to_vcf.out.vcf, bcf_to_vcf.out.index, assembly, run)
//                quality_metrics = QC1.concat(QC2, annotation_table_merged.out.vep_merged_stat).collect()
//                multiqc_pop(quality_metrics, assembly, run, SNV)

// Step reoved and hopefully to add to Hail
//              gnomad_frequency_table(gnomad_SNV_vcf, gnomad_SNV_index, chr)
//              split_vcf_by_chr(bcf_to_vcf.out.vcf, assembly, batch, run, chr)
//                SNV_data_organization(gnomad_frequency_table.out.collect(), split_vcf_by_chr.out.vcf_onechr, annotation_table_merged.out.annot_table_merged_R.collect(), assembly, run, sample_QC.out.sample_QC_file)

