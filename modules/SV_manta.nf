// Nextflow process
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// SV calling. Call SV using manta
// Change the header of the vcf and index the compressed vcf


process SV_manta {
	tag "${bam.simpleName}"
    	
	publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/SV/Sample/manta", mode: 'copyNoFollow'
	
	input:
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
	if [ -a $params.outdir_ind/${assembly}/*/${run}/SV/Sample/manta/${bam.simpleName}_diploidSV.vcf.gz ]; then
		manta_vcf=\$(find $params.outdir_ind/${assembly}/*/${run}/SV/Sample/manta/ -name ${bam.simpleName}_diploidSV.vcf.gz)
		manta_index=\$(find $params.outdir_ind/${assembly}/*/${run}/SV/Sample/manta/ -name ${bam.simpleName}_diploidSV.vcf.gz.tbi)
		ln -s \$manta_vcf .
		ln -s \$manta_index .
	else
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
	fi
	"""
}
