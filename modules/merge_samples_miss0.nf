// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// Merge all the samples (without joint calling) and create a vcf files with all the calls for each participant
// Set the missing varaints to 0 as no joint calling is done for the MEI
// Index the vcf file

process merge_samples_miss0 {
	label 'conda_annotate'

        publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/${var_type}/", mode: 'copy'

        input :
        file list_vcf
	val assembly
	val batch
	val run
	val var_type

        output:
        path '*.vcf.gz', emit : vcf

        script :
        """
        bcftools merge --missing-to-ref -m none -l ${list_vcf} -o ${var_type}_${run}.vcf
	bcftools view -i "MAC >=1" ${var_type}_${run}.vcf -O z -o ${var_type}_${run}.vcf.gz
        """
}
