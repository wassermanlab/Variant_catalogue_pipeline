// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// This step filters out blacklisted sites containing unwanted artifacts.


process MT_FilterOut_sites {
        tag "${MT_trimmed.simpleName}"

        publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/MT/Sample/", mode: 'copy'

        input :
        file ref_genome_MT
        file ref_genome_MT_index
	file ref_genome_MT_dict
	file (MT_trimmed)
        file MT_trimmed_index
	file blacklist_sites_hg38_MT_file
	file blacklist_sites_hg38_MT_file_index
	val assembly
	val batch
	val run

        output :
        path '*_filtered_sites.vcf.gz', emit : vcf
	path '*_filtered_sites.vcf.gz.tbi', emit : index
	path '*_MT_Step2_participant_data.tsv', emit : sample_MT_Step2_participant_data
	path '*_list.txt', emit : Sample_list

        script :
        """
        sample_name=\$(echo ${MT_trimmed.simpleName} | sed 's/_.*//' )
	
	gatk VariantFiltration \
	-R Homo_sapiens_assembly38.chrM.fasta \
	-V ${MT_trimmed.simpleName}.vcf.gz \
	-O ${MT_trimmed.simpleName}_filtered_sites.vcf.gz \
	--mask-name "GATK_artifact" \
	--mask ${blacklist_sites_hg38_MT_file}

	echo "\${sample_name}\t\${sample_name}\t$params.outdir_ind/${assembly}/${batch}/${run}/MT/Sample/${MT_trimmed.simpleName}_filtered_sites.vcf.gz" > \${sample_name}_MT_Step2_participant_data.tsv
	echo "\${sample_name}" > \${sample_name}_list.txt
        """
}

