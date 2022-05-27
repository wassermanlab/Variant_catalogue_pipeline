// Nextflow process
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// SV calling
// STEP TO DESCRIBE

process SV_paragraph_duphold {
        tag "${bam.simpleName}"

//	Issue with the script from paragraph   [E::idx_find_and_load] Could not retrieve index file for 'paragraph_output/variants.vcf.gz'
// 	According to github, workaround is to modify multigrmpy.py : https://github.com/Illumina/paragraph/issues/59
	errorStrategy 'ignore' // TODO: change after debugging is done

        publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/SV/Sample/paragraph", mode: 'copy', pattern : '*_genotypes_setid.vcf.gz*'

	input:
	tuple(path(site_vcf), path(site_vcf_index))
	file bam
	file bai
	file reference
	file reference_index
	val assembly
	val batch
	val run

	output:
        path '*_genotypes_setid.vcf.gz.tbi', emit : index 
	path '*_genotypes_setid.vcf.gz', emit : vcf	

	script:
	output_file = "${bam.simpleName}.paragraph.vcf.gz"
	"""
	dp=\$(tiwih meandepth $bam)
	tsample=\$(tiwih samplename $bam)
        echo "\$tsample" > sample.txt
	echo "id\tpath\tdepth\tread length" > sample.manifest
	echo "\$tsample\t$bam\t\$dp\t150" >> sample.manifest
	M=\$((dp * 5))
	cat sample.manifest
    
	# this is the main paragraph entrypoint
	multigrmpy.py -i $site_vcf \
        -m sample.manifest \
        -r $reference \
        -o paragraph_output \
        -t ${task.cpus} \
        -M \$M

	ANNOTATEVARIANTS_INSTALL=/mnt/common/WASSERMAN_SOFTWARE/AnnotateVariants/
	source \$ANNOTATEVARIANTS_INSTALL/opt/miniconda3/etc/profile.d/conda.sh
	conda activate \$ANNOTATEVARIANTS_INSTALL/opt/AnnotateVariantsEnvironment
	
	bcftools index -t paragraph_output/genotypes.vcf.gz

	bcftools view -S sample.txt -O z -o ${bam.simpleName}_genotypes.vcf.gz paragraph_output/genotypes.vcf.gz
        bcftools annotate --set-id '%CHROM\\_%POS\\_%SVTYPE\\_%SVLEN' -O z -o ${bam.simpleName}_genotypes_setid.vcf.gz ${bam.simpleName}_genotypes.vcf.gz

	bcftools index -t ${bam.simpleName}_genotypes_setid.vcf.gz

	# duphold adds depth annotations looking at coverage fold-change around Svs
	#duphold -d -v t/genotypes.vcf.gz -b $bam -f $reference -t 4 -o $output_file
	#bcftools index --threads 3 $output_file
	"""
}
