// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// MT. Identify Duplicate Reads using MarkDuplicates
// This step identifies and tags duplicate reads in the aligned BAM files.



process MarkDuplicates {
        tag "${bam_MT.baseName}"

        input :
        file bam_MT
	file bai_MT

        output :
        path '*marked_duplicates.bam', emit : bam
        path '*marked_duplicates.bam.bai', emit : bai

        script:
        """
        singularity exec -B /mnt/scratch/SILENT/Act3/ -B /mnt/common/SILENT/Act3/ /mnt/common/SILENT/Act3/singularity/gatk4-4.2.0.sif \
	gatk MarkDuplicates \
	I=${bam_MT.baseName}.bam \
	O=${bam_MT.baseName}_marked_duplicates.bam \
	M=${bam_MT.baseName}_marked_duplicates_metrics.txt

	ANNOTATEVARIANTS_INSTALL=/mnt/common/WASSERMAN_SOFTWARE/AnnotateVariants/
        source \$ANNOTATEVARIANTS_INSTALL/opt/miniconda3/etc/profile.d/conda.sh
        conda activate \$ANNOTATEVARIANTS_INSTALL/opt/AnnotateVariantsEnvironment

	samtools index ${bam_MT.baseName}_marked_duplicates.bam
        """
}
