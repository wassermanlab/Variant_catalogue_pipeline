// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// Merge all the samples (without joint calling) and create a vcf files with all the calls for each participant
// Index the vcf file

process merge_samples {

        publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/${var_type}/", mode: 'copy'

        input :
        file list_vcf
	val assembly
	val batch
	val run
	val var_type

        output:
        path '*.vcf.gz', emit : vcf
	path '*.vcf.gz.tbi', emit : index

        script :
        """
	ANNOTATEVARIANTS_INSTALL=/mnt/common/WASSERMAN_SOFTWARE/AnnotateVariants/
	source \$ANNOTATEVARIANTS_INSTALL/opt/miniconda3/etc/profile.d/conda.sh
	conda activate \$ANNOTATEVARIANTS_INSTALL/opt/AnnotateVariantsEnvironment

        bcftools merge -m none -l ${list_vcf} -o ${var_type}_${run}.vcf
	bcftools view -i "MAC >=1" ${var_type}_${run}.vcf -O z -o ${var_type}_${run}.vcf.gz

        singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/gatk4-4.2.0.sif \
        gatk --java-options "-Xmx4G" \
	IndexFeatureFile \
        -I ${var_type}_${run}.vcf.gz
        """
}
