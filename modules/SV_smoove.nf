// Nextflow process
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// SV calling. Call variants using smoove
// compress, change the header and index the compressed file


process SV_smoove {
	tag "${bam.simpleName}"

	publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/SV/Sample/smoove", mode: 'copy'
	
	input:
	file bam
	file bai
	file reference
	file reference_index
	val assembly
	val batch
	val run

	output:
	tuple(file("${bam.simpleName}_smoove.vcf.gz"),file("${bam.simpleName}_smoove.vcf.gz.tbi"),val(bam.simpleName))	
	
	script:
	"""
	echo ${bam.simpleName}
	sample_name=\$(echo ${bam.simpleName} | cut -d _ -f 1)
	echo \$sample_name > sample.txt

	smoove call \
	--outdir . \
	--name ${bam.simpleName} \
	--fasta ${reference}\
	${bam}
  	##-p ${task.cpus} 
	
	bcftools view -O u -o ${bam.simpleName}.R.bcf ${bam.simpleName}-smoove.vcf.gz
	bcftools sort --temp-dir $params.outdir_ind/${assembly}/${batch}/${run}/SV/TMP  -m 2G -O z -o ${bam.simpleName}-smoove.vcf.gz  ${bam.simpleName}.R.bcf
	bcftools reheader -s sample.txt ${bam.simpleName}-smoove.vcf.gz > ${bam.simpleName}_smoove.vcf.gz
	bcftools index --tbi ${bam.simpleName}_smoove.vcf.gz
	"""
}
