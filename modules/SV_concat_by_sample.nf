// Nextflow process
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// SV calling. Concatenante the SVs called by different callers for each sample



process SV_concat_by_sample {
	tag "${sample_name}"

	publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/SV/Sample/Concat_by_sample"

	input:
        tuple(path(vcfs), path(indexes), val(sample_name))
	val assembly
	val batch
	val run	

	output:
	file("${output_file}")

	script:
	output_file = "${sample_name}.concat-svs.vcf"
	"""
	if [ -a $params.outdir_ind/${assembly}/*/${run}/SV/Sample/Concat_by_sample/${sample_name}.concat-svs.vcf ]; then
		concat_vcf=\$(find $params.outdir_ind/${assembly}/*/${run}/SV/Sample/Concat_by_sample/  -name ${sample_name}.concat-svs.vcf) 
		ln -s \$concat_vcf .
	else
		bcftools concat -a -O v -o ${output_file} *.vcf.gz
	fi
	"""
}
