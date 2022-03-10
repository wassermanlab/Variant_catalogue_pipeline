// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// Merge stats files for output VCFs
// The output file is necessary to run GATK FilterMutectCalls 

process MT_Merge_stat_file {
        tag "${MT_call_variants_stat.simpleName}"

        input :
        file MT_call_variants_stat
	file MT_call_variants_shifted_stat

        output :
        file '*'

        script :
        """
        echo ${MT_call_variants_stat.simpleName}
        sample_name=\$(echo ${MT_call_variants_stat.simpleName} | cut -d _ -f 1)
        echo \$sample_name
        
        gatk MergeMutectStats \
        -stats ${MT_call_variants_stat} \
        -stats \${sample_name}_sorted_chrM_Homo_sapiens_assembly38.chrM.shifted_by_8000_bases_marked_duplicates_Mutect2.vcf.gz.stats \
        -O \${sample_name}_MT_merged.stats
	"""
}
