// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// Create a file with the list of variants and their annotation
// Same process for SNV and MT variants (using vep, obtaining a tsv)
// The last line remove the # in front of "Uploaded_varaintion", which is necessary for downstream analysis

process annotation_table_merged {
	label 'conda_annotate'
	tag "${chr}"

    publishDir "$params.outdir_pop/${assembly}/${run}/${var_type}/VEP_annotation/", mode: 'copy', pattern : '*_annotation_table_merged*'
	publishDir "$params.outdir_pop/${assembly}/${run}/QC/${var_type}/", mode: 'copy', pattern : '*_VEP_merged_stats*'

        input :
        file vcf
		file index
		val vep_cache_merged
		val vep_cache_merged_version
		val assembly
		val run
		val assembly_VEP
		file CADD_1_6_whole_genome_SNVs
		file CADD_1_6_whole_genome_SNVs_index
		file CADD_1_6_InDels
		file CADD_1_6_InDels_index
		file spliceai_snv
		file spliceai_snv_index
		file spliceai_indel
        file spliceai_indel_index
		each chr
		val var_type 
		file reference
		file dir_plugin

        output :

        path '*_VEP_merged_stats*', emit : vep_merged_stat
		path '*.vcf', emit :  annotation_vcf

        script :
        """
	vep \
        -i ${vcf} \
        -o ${vcf.simpleName}_${var_type}_annotation_table_merged_${chr}.vcf \
	--vcf \
	--chr ${chr}  \
	--offline \
	--merged \
        --assembly $assembly_VEP \
        --cache \
        --dir_cache ${vep_cache_merged} \
	--cache_version ${vep_cache_merged_version} \
	--use_transcript_ref \
	--fasta $reference \
        --distance 0 \
	--symbol \
	--biotype \
	--transcript_version \
	--variant_class \
	--polyphen b \
	--sift b \
	--hgvs \
        --check_existing \
        --var_synonyms \
	--tsl \
        --dir_plugins ${dir_plugin} \
	--plugin CADD,$CADD_1_6_whole_genome_SNVs,$CADD_1_6_InDels \
        --plugin SpliceAI,snv=${spliceai_snv},indel=${spliceai_indel} \
	--stats_file ${vcf.simpleName}_${chr}_VEP_merged_stats
	"""
}
