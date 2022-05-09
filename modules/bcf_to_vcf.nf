// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// SNV Calling. 
// Split the multiallelic varaints (norm step) and transform the bcf into a vcf 
// Rename the varaints and compress the vcf into a vcf.gz
// Index the compressed vcf

process bcf_to_vcf {
        publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/SNV/", mode: 'copy'
//        publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/vcf_pre_hail/", mode: 'copy'

	input :
	file bcf_file
	val assembly
	val batch
	val run

	output :
	path '*.vcf.gz', emit : vcf
	path '*.vcf.gz.tbi', emit : index
	

	script :
	"""
	ANNOTATEVARIANTS_INSTALL=/mnt/common/WASSERMAN_SOFTWARE/AnnotateVariants/
	source \$ANNOTATEVARIANTS_INSTALL/opt/miniconda3/etc/profile.d/conda.sh
	conda activate \$ANNOTATEVARIANTS_INSTALL/opt/AnnotateVariantsEnvironment

##	bcftools view ${bcf_file} | bgzip -c > ${bcf_file.simpleName}.vcf.gz
	bcftools norm -m -any -o ${bcf_file.simpleName}_norm.vcf ${bcf_file}
	bcftools annotate --set-id '%CHROM\\_%POS\\_%REF\\_%FIRST_ALT' -O z -o ${bcf_file.simpleName}.vcf.gz ${bcf_file.simpleName}_norm.vcf
	
        singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/gatk4-4.2.0.sif \
	gatk --java-options "-Xmx4G" \
	IndexFeatureFile \
        -I ${bcf_file.simpleName}.vcf.gz	
	"""
}
