// Nextflow sub-workflow
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the sub-workflow goal and characteristics :
// Call the complex variants : Structural varaints (SV), Mobile Element Insertions (MEI) and Short Tendem Repeats (STR)


// Load the modules for the SV workflow

include { SV_smoove } from "./../modules/SV_smoove"
include { SV_manta } from "./../modules/SV_manta"
include { SV_concat_by_sample } from "./../modules/SV_concat_by_sample"
include { SV_jasmine } from "./../modules/SV_jasmine"
include { SV_paragraph_duphold } from "./../modules/SV_paragraph_duphold"

include { samtools_fixmate } from "./../modules/samtools_fixmate"
include { expansion_hunter } from "./../modules/expansion_hunter"
include { melt } from "./../modules/melt"

include { list_vcfs_txt as SV_vcfs_txt; list_vcfs_txt as STR_vcfs_txt; list_vcfs_txt as MEI_vcfs_txt } from "./../modules/list_vcfs_txt"
include { merge_samples as SV_merge_samples} from "./../modules/merge_samples"
include { merge_samples_miss0 as MEI_merge_samples; merge_samples_miss0 as STR_merge_samples } from "./../modules/merge_samples_miss0"
include { annotation_table_merged as SV_annotation; annotation_table_merged as MEI_annotation } from "./../modules/annotation_table_merged"

include { SV_data_organization } from "./../modules/SV_data_organization"
include { MEI_data_organization } from "./../modules/MEI_data_organization"
include { STR_data_organization } from "./../modules/STR_data_organization"

include { Hail_MEI_QC } from "./../modules/Hail_MEI_QC"
include { Hail_STR } from "./../modules/Hail_STR"
include { Hail_SV_QC } from "./../modules/Hail_SV_QC"

// SV workflow

workflow SV {

	// Load the parameters and files
	run             			= params.run
	batch           			= params.batch
	assembly        			= params.assembly
	reference       			= file (params.ref)
	reference_index 			= file (params.ref_index)
        transposon_file				= file (params.transposon_file)
        cr_bed 					= file (params.cr_bed)
        cr_bed_index                            = file (params.cr_bed_index)
	genes_file				= file (params.genes_file)
	variant_catalog				= file (params.variant_catalog)
	chr					= params.chrom
	SV					= params.SV
	STR					= params.STR
	MEI					= params.MEI
//	ExpansionHunterDenovo			= file (params.ExpansionHunterDenovo)
//	manifest_STR				= file (params.manifest_STR)

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
	dir_plugin				= file (params.dir_plugin)
        severity_table                          = file (params.severity_table)

	// Workflow start here
	take : 
		bam
		bai
		sample_sex_file

	main :
		//Structural Varaints (SV)
                // Sample specific (Do not need to be run for a previously processed sample)
		sm = SV_smoove(bam, bai, reference, reference_index, assembly, batch, run)
		mr = SV_manta(bam, bai, reference, reference_index, cr_bed, cr_bed_index, assembly, batch, run)
		sv_groups = mr.concat(sm) | groupTuple(by: 2)
		svs = SV_concat_by_sample(sv_groups, assembly, batch, run) | collect
		sv_merged = SV_jasmine(svs, reference, reference_index, assembly, batch, run)
		genotyped = SV_paragraph_duphold(sv_merged, bam, bai, reference, reference_index, assembly, batch, run)

                // Aggregated steps (Need to be run everytime a new sample is added to the cohort)
		SV_vcfs_txt(SV_paragraph_duphold.out.vcf.collect(), assembly, batch, run, SV)
		SV_merge_samples(SV_vcfs_txt.out, assembly, batch, run, SV)
                Hail_SV_QC (SV_merge_samples.out.vcf, sample_sex_file, assembly, batch, run)
		SV_annotation(Hail_SV_QC.out.vcf_SV_filtered_frequ_only, Hail_SV_QC.out.index_SV_filtered_frequ_only, vep_cache_merged, vep_cache_merged_version, assembly, run, assembly, CADD_1_6_whole_genome_SNVs, CADD_1_6_whole_genome_SNVs_index, CADD_1_6_InDels, CADD_1_6_InDels_index, spliceai_snv, spliceai_snv_index, spliceai_indel, spliceai_indel_index, chr, SV, reference, dir_plugin)

                SV_data_organization(SV_annotation.out.annotation_vcf, assembly, run, SV, severity_table)


		//Short Tandem Repeats (STR)
                // Sample specific (Do not need to be run for a previously processed sample)
		expansion_hunter(bam, bai, reference, reference_index, variant_catalog, assembly, batch, run)

                // Aggregated steps (Need to be run everytime a new sample is added to the cohort)
		STR_vcfs_txt(expansion_hunter.out.vcf.collect(), assembly, batch, run, STR)
  		STR_merge_samples(STR_vcfs_txt.out, assembly, batch, run, STR)
		Hail_STR (STR_merge_samples.out.vcf, sample_sex_file, assembly, batch, run) 
                STR_data_organization(STR_merge_samples.out.vcf, variant_catalog, assembly, run, STR)



		// Mobile Element Insertions (MEIs)
                // Sample specific (Do not need to be run for a previously processed sample)
		samtools_fixmate(bam, bai, assembly, batch, run)
		melt(samtools_fixmate.out.samples_fixmate_bam, samtools_fixmate.out.samples_fixmate_bam_index, reference, reference_index, transposon_file, genes_file, assembly, batch, run)
                
		// Aggregated steps (Need to be run everytime a new sample is added to the cohort)
		MEI_vcfs_txt(melt.out.vcf.collect(), assembly, batch, run, MEI)
		MEI_merge_samples(MEI_vcfs_txt.out, assembly, batch, run, MEI)
                Hail_MEI_QC (MEI_merge_samples.out.vcf, sample_sex_file, assembly, batch, run)
                MEI_annotation(Hail_MEI_QC.out.vcf_MEI_filtered_frequ_only, Hail_MEI_QC.out.index_MEI_filtered_frequ_only, vep_cache_merged, vep_cache_merged_version, assembly, run, assembly, CADD_1_6_whole_genome_SNVs, CADD_1_6_whole_genome_SNVs_index, CADD_1_6_InDels, CADD_1_6_InDels_index, spliceai_snv, spliceai_snv_index, spliceai_indel, spliceai_indel_index, chr, MEI, reference, dir_plugin)

                 MEI_data_organization(MEI_annotation.out.annotation_vcf, assembly, run, MEI)
}
