// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// SNV Caling. GLnexus to do joint variant calling
// GLnexus also include a varaint filtering step, that according to publications, is as good as GATK VQSR, so no additional step is needed.

process GLnexus_cli {

	input :
	file list_gvcf
	val run
        
	output :
	path '*.bcf'

	script :
	"""
	glnexus_cli \
	--config DeepVariantWGS \
	--mem-gbytes 128 \
        --threads ${task.cpus} \
	--list ${list_gvcf} > DeepVariant_GLnexus_${run}.bcf
	"""
}

