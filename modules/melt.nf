// Nextflow process
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// SV calling. Call Mobile Element Insertions (MEI) using MELT
// Rename varaints, compress the vcf and index the compressed vcf

process melt {
	label 'conda_annotate'
	tag "${bam.simpleName}"

	publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/MEI/Sample/"//, mode: 'copyNoFollow'
		
	input:
	file bam
	file bai
	file reference
	file reference_index
	file transposon_file
	file genes_file
	val assembly
	val batch
	val run

	output:
	path '*_mei.vcf.gz', emit : vcf	
	path '*_mei.vcf.gz.tbi', emit : vcf_index

	script:
	"""
		sample_name=\$(echo ${bam.simpleName} | cut -d _ -f 1)
		# Resort the bam file by query name for samtools fixmate (coordiante-sorted bam does not work)
		samtools sort -n -O BAM -@ 20  ${bam} > ${bam.SimpleName}_nsorted.bam

		# Samtools fixmate will add MQ tags
		samtools fixmate -m -O BAM -@ 20  ${bam.SimpleName}_nsorted.bam  ${bam.SimpleName}_fixmate.bam

		# Now sort bam file by coordinates to resume the pipeline 
		samtools sort  -@ 20  ${bam.SimpleName}_fixmate.bam -o ${bam.SimpleName}_fixmate_ordered.bam

		# index the sorted bam file
		samtools index  -@ 20  ${bam.SimpleName}_fixmate_ordered.bam
		rm *_nsorted.bam* *_fixmate.bam

        sample_name=\$(echo ${bam.simpleName} | cut -d _ -f 1)
		mkdir -p \${sample_name}
	
		java -Xmx8G -jar ${params.Melt_dir}/MELT.jar Single \
		-b hs37d5/NC_007605 \
		-t ${transposon_file}  \
		-h ${reference} \
		-bamfile ${bam.SimpleName}_fixmate_ordered.bam \
		-w \${sample_name} \
		-n ${genes_file}

		# fix issues with  MELT vcfs
		for name in {ALU,SVA,LINE1};do bcftools annotate -x FMT/GL \${sample_name}/\${name}.final_comp.vcf > \${name}.vcf;done
		for name in {ALU,SVA,LINE1};do bcftools view -O u -o \${name}.bcf \${name}.vcf;done
		for name in {ALU,SVA,LINE1};do bcftools sort  -m 2G -O z -o \${name}.vcf.gz \${name}.bcf;done
		for name in {ALU,SVA,LINE1};do bcftools index --tbi \${name}.vcf.gz;done
		bcftools concat -a -Oz  -o \${sample_name}_mei_noID.vcf.gz *vcf.gz
		bcftools annotate --set-id '%CHROM\\_%POS\\_%SVTYPE\\_%SVLEN' -O z -o \${sample_name}_mei.vcf.gz \${sample_name}_mei_noID.vcf.gz
		bcftools index --tbi \${sample_name}_mei.vcf.gz
		rm *fixmate_ordered.bam*
	"""
}

