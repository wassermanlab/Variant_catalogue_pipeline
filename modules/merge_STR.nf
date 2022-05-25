// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// Merge all the samples (without joint calling) and create a vcf files with all the calls for each participant
// Index the vcf file

process merge_STR {

        publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/STR/", mode: 'copy'

        input :
//        file manifest
	file reference
	val assembly
	val batch
	val run

        output:
        path '*.json', emit : output_json

        script :
        """
	/mnt/common/Precision/ExpansionHunterDenovo/ExpansionHunterDenovo-v0.9.0-linux_x86_64/bin/ExpansionHunterDenovo merge\
	--reference ${reference} \
	--manifest /mnt/scratch/SILENT/Act3/Processed/Workflow/Solenne/IBVL_pipeline/manifest.txt \
	--output-prefix ${run}
        """
}
