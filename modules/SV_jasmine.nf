// Nextflow process
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// SV calling. 
// DESCRIPTION OF THIS STEP TO DO



process SV_jasmine {
	publishDir "$params.outdir_ind/${assembly}/${batch}/${run}/SV/jasmine", mode: 'copy'

	input:
        val(sample_vcfs)
	file reference
	file reference_index
	val assembly
	val batch
	val run

	output: 
	tuple(file("${output_file}"), file("${output_file}.tbi"))

	script:
	output_file = "jasmine.merged.vcf.gz"
	tmp_dir = "$params.outdir_ind/${assembly}/${batch}/${run}/SV"
	file("$workDir/vcfs.list").withWriter { fh ->
		sample_vcfs.each { vcf ->
                	// write .vcf even though it's really .vcf.gz since jasmine doesn't accept gz
                	// and we change the file below.
                	fh.write(vcf.toString()); fh.write("\n")
		}
	}

	// we don't merge if we only have a single sample.
	if(sample_vcfs.size() > 1) {	
		"""
        	# jasmine can't do gzip.
        	jasmine  -Xmx6g \
		--threads ${task.cpus} \
		--output_genotypes \
		--allow_intrasample \
		--file_list=${workDir}/vcfs.list out_file=${output_file} \
		genome_file=$reference

		#Removed option --dup_to_ins in jasmine because it caused error for the next step (paraGRAPH) , to consider adding again in future version.
#                --dup_to_ins  \

        
		# NOTE: removing BNDs and setting min start to > 150 as paragraph fails if start < readlength
        	tiwih setsvalt --drop-bnds --inv-2-ins -o ${output_file}.tmp.vcf.gz $reference $output_file
        	bcftools sort --temp-dir ${tmp_dir}  -m 2G -O z -o ${output_file} ${output_file}.tmp.vcf.gz
        	bcftools index --tbi $output_file
        	"""
	} else {
        	"""
        	tiwih setsvalt --drop-bnds -o ${output_file} $reference ${sample_vcfs[0]}
        	tabix $output_file
        	"""
	}

}

