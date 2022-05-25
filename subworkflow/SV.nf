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
include { merge_samples as SV_merge_samples; merge_samples as STR_merge_samples; merge_samples as MEI_merge_samples } from "./../modules/merge_samples"

include { merge_STR } from "./../modules/merge_STR"
include { SV_split_vcf_by_chr as SV_split_vcf_by_chr; SV_split_vcf_by_chr as MEI_split_vcf_by_chr } from "./../modules/SV_split_vcf_by_chr"
include { annotation_table_merged as SV_annotation; annotation_table_merged as MEI_annotation } from "./../modules/annotation_table_merged"
include { SV_data_organization } from "./../modules/SV_data_organization"
include { MEI_data_organization } from "./../modules/MEI_data_organization"
include { STR_data_organization } from "./../modules/STR_data_organization"


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

        vep_cache_merged                        = file (params.vep_cache_merged)
        vep_cache_merged_version                = params.vep_cache_merged_version
        CADD_1_6_whole_genome_SNVs              = file (params.CADD_1_6_whole_genome_SNVs)
        CADD_1_6_whole_genome_SNVs_index        = file (params.CADD_1_6_whole_genome_SNVs_index)
        CADD_1_6_InDels                         = file (params.CADD_1_6_InDels)
        CADD_1_6_InDels_index                   = file (params.CADD_1_6_InDels_index)

	// Workflow start here
	take : 
		bam
		bai
		sample_sex_file

	main :
		sm = SV_smoove(bam, bai, reference, reference_index, assembly, batch, run)
		mr = SV_manta(bam, bai, reference, reference_index, cr_bed, cr_bed_index, assembly, batch, run)
		sv_groups = mr.concat(sm) | groupTuple(by: 2)
		svs = SV_concat_by_sample(sv_groups, assembly, batch, run) | collect
		sv_merged = SV_jasmine(svs, reference, reference_index, assembly, batch, run)
		genotyped = SV_paragraph_duphold(sv_merged, bam, bai, reference, reference_index, assembly, batch, run)
		SV_vcfs_txt(SV_paragraph_duphold.out.vcf.collect(), assembly, batch, run, SV)
		SV_merge_samples(SV_vcfs_txt.out, assembly, batch, run, SV)

		expansion_hunter(bam, bai, reference, reference_index, variant_catalog, assembly, batch, run)
		STR_vcfs_txt(expansion_hunter.out.vcf.collect(), assembly, batch, run, STR)
  		STR_merge_samples(STR_vcfs_txt.out, assembly, batch, run, STR)

		samtools_fixmate(bam, bai, assembly, batch, run)
		melt(samtools_fixmate.out.samples_fixmate_bam, samtools_fixmate.out.samples_fixmate_bam_index, reference, reference_index, transposon_file, genes_file, assembly, batch, run)
		MEI_vcfs_txt(melt.out.vcf.collect(), assembly, batch, run, MEI)
		MEI_merge_samples(MEI_vcfs_txt.out, assembly, batch, run, MEI)

                MEI_split_vcf_by_chr(MEI_merge_samples.out.vcf, assembly, batch, run, chr, MEI)
                MEI_annotation(MEI_merge_samples.out.vcf, MEI_merge_samples.out.index, vep_cache_merged, vep_cache_merged_version, assembly, run, assembly, CADD_1_6_whole_genome_SNVs, CADD_1_6_whole_genome_SNVs_index, CADD_1_6_InDels, CADD_1_6_InDels_index, chr, MEI)
                MEI_data_organization(MEI_split_vcf_by_chr.out.vcf_onechr, MEI_annotation.out.annot_table_merged_R.collect(), assembly, run, MEI, sample_sex_file)
                STR_data_organization(STR_merge_samples.out.vcf, variant_catalog, assembly, run, STR)
                SV_split_vcf_by_chr(SV_merge_samples.out.vcf, assembly, batch, run, chr, SV)
                SV_annotation(SV_merge_samples.out.vcf, SV_merge_samples.out.index, vep_cache_merged, vep_cache_merged_version, assembly, run, assembly, CADD_1_6_whole_genome_SNVs, CADD_1_6_whole_genome_SNVs_index, CADD_1_6_InDels, CADD_1_6_InDels_index, chr, SV)
                SV_data_organization(SV_split_vcf_by_chr.out.vcf_onechr, SV_annotation.out.annot_table_merged_R.collect(), assembly, run, SV, sample_sex_file)
}
