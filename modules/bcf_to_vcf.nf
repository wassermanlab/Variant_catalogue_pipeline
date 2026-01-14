// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developed to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// SNV Calling. 
// Optionally subset VCF by specified populations
// Split the multiallelic variants (norm step) and transform the bcf into a vcf 
// Rename the variants and compress the vcf into a vcf.gz
// Index the compressed vcf

process bcf_to_vcf {
	label 'conda_annotate'
        publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/SNV/", mode: 'copy'

	input :
	file bcf_file
        val assembly
        val batch
        val run
        file ref
        file ref_index
        path sample_assignments
        path pop_list

        output:
        path '*_norm.vcf.gz', emit : vcf                 
        path '*GLnexus_output.vcf.gz'
        path '*.vcf.gz.tbi'

	script :
	"""
        # subset the VCF based on the populations of interest
        grep -f ${pop_list} ${sample_assignments} > subset_assignments.txt
        cut -d',' -f1 subset_assignments.txt > sample_subset_list.txt

        # subset the population vcf
        bcftools view -S sample_subset_list.txt ${bcf_file} -Oz -o ${bcf_file.simpleName}_GLnexus_output.vcf.gz

        # index
        bcftools index -t ${bcf_file.simpleName}_GLnexus_output.vcf.gz

        # normalize/left align and split multi-allelic variants
        bcftools norm -m -any -Oz -o ${bcf_file.simpleName}_norm_int.vcf.gz \
            -f ${ref} ${bcf_file.simpleName}_GLnexus_output.vcf.gz
        bcftools index -t  ${bcf_file.simpleName}_norm_int.vcf.gz
        bcftools annotate --set-id '%CHROM\\_%POS\\_%REF\\_%FIRST_ALT' -O z -o ${bcf_file.simpleName}_norm.vcf.gz \
            ${bcf_file.simpleName}_norm_int.vcf.gz

	"""
}
