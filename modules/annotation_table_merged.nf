// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// Create a file with the list of variants and their annotation
// Same process for SNV and MT variants (using vep, obtaining a tsv)
// The last line remove the # in front of "Uploaded_varaintion", which is necessary for downstream analysis

process annotation_table_merged {
	tag "${chr}"

        publishDir "$params.outdir_pop/${assembly}/${run}/VEP_annotation/", mode: 'copy', pattern : '*_annotation_table_merged.tsv'
	publishDir "$params.outdir_pop/${assembly}/${run}/QC/${vcf.simpleName}/", mode: 'copy', pattern : '*_VEP_stats*'

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
	each chr

        output :
        path '*_annotation_table_merged.tsv', emit :  annotation_table_merged
        path '*_VEP_merged_stats*', emit : vep_merged_stat
	path '*_annotation_table_merged_nohash.tsv', emit : annot_table_merged_R

        script :
        """
	ANNOTATEVARIANTS_INSTALL=/mnt/common/WASSERMAN_SOFTWARE/AnnotateVariants/
        source \$ANNOTATEVARIANTS_INSTALL/opt/miniconda3/etc/profile.d/conda.sh
	conda activate \$ANNOTATEVARIANTS_INSTALL/opt/AnnotateVariantsEnvironment

	vep \
        -i ${vcf} \
        -o ${vcf.simpleName}_${chr}_annotation_table_merged.tsv \
	--chr ${chr}  \
	--offline \
	--merged \
        --assembly $assembly_VEP \
        --cache \
        --dir_cache ${vep_cache_merged} \
	--cache_version ${vep_cache_merged_version} \
	--use_transcript_ref \
	--fasta /mnt/common/DATABASES/REFERENCES/GRCh37/GENOME/GRCh37-lite.fa \
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
        --tab \
        --dir_plugins /mnt/common/SILENT/Act3/VEP/Plugins/ \
        --plugin CADD,$CADD_1_6_whole_genome_SNVs,$CADD_1_6_InDels \
        --stats_file ${vcf.simpleName}_${chr}_VEP_merged_stats

	sed 's/#Uploaded_variation/Uploaded_variation/g' ${vcf.simpleName}_${chr}_annotation_table_merged.tsv > ${vcf.simpleName}_${chr}_annotation_table_merged_nohash.tsv
        """
}
