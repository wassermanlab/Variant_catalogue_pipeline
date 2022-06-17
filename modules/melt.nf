// Nextflow process
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// SV calling. Call Mobile Element Insertions (MEI) using MELT
// Rename varaints, compress the vcf and index the compressed vcf

process melt {
	tag "${bam.simpleName}"

	publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/MEI/Sample/", mode: 'copyNoFollow'
	
//	module = "BCCHR/Java/1.8.0_231"
//	conda  "/home/BCRICWH.LAN/Mohammed.Mohammed/miniconda3/envs/sv"   
	
	input:
//	tuple(val(sample_name), path(bam), path(index))
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

	if [ -a $params.outdir_ind/${assembly}/*/${run}/MEI/Sample/\${sample_name}_mei.vcf.gz ]; then
		melt_vcf=\$(find $params.outdir_ind/${assembly}/*/${run}/MEI/Sample/ -name \${sample_name}_mei.vcf.gz)
		melt_index=\$(find $params.outdir_ind/${assembly}/*/${run}/MEI/Sample/ -name \${sample_name}_mei.vcf.gz.tbi)
		ln -s \$melt_vcf .
		ln -s \$melt_index .
	else
        	# Unload bcchr, and load cvmfs
        	# unload_bcchr
        	source /cm/shared/BCCHR-apps/env_vars/unset_BCM.sh
        	# load cvmfs
        	source /cvmfs/soft.computecanada.ca/config/profile/bash.sh

        	module load StdEnv/2020
        	module load vcftools
        	module load bcftools
		module load bowtie2

		mkdir -p \${sample_name}
	
		java -Xmx8G -jar ${params.Melt_dir}/MELT.jar Single \
		-b hs37d5/NC_007605 \
		-t ${transposon_file}  \
		-h ${reference} \
		-bamfile $bam \
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
	fi
	"""
}

