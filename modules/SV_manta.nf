process SV_manta {
	tag "${bam.simpleName}"
    	
	publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/SV/manta", mode: 'copy'
	
	input:
//	tuple(val(sample_name), path(bam), path(index))
	file bam
	file bai
	file reference
	file reference_index
	file cr_bed
	file cr_bed_index
	val assembly
	val batch
	val run
	
	output:
	tuple(file("${bam.simpleName}_diploidSV.vcf.gz"),file("${bam.simpleName}_diploidSV.vcf.gz.tbi"), val(bam.simpleName))
   
	script:
	"""
        echo ${bam.simpleName}
        sample_name=\$(echo ${bam.simpleName} | cut -d _ -f 1)
        echo \$sample_name > sample.txt

	configManta.py \
	--bam ${bam} \
	--referenceFasta ${reference} \
	--runDir . \
	--callRegions ${cr_bed}

	python2 ./runWorkflow.py \
	-j ${task.cpus} \
	-m local
	
	bcftools reheader -s sample.txt results/variants/diploidSV.vcf.gz > ${bam.simpleName}_diploidSV.vcf.gz
	bcftools index -f --tbi ${bam.simpleName}_diploidSV.vcf.gz
	"""
}
