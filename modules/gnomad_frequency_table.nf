// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// Extract the useful information from the gnomAD table to make it smaller and load it in R
// Followed by a step splitting the output by chr?
// To specify a chr -L ${chr}

process gnomad_frequency_table {
	tag "${gnomad_vcf}"

	publishDir "$params.gnomad_SNV_frequency_outdir", mode: 'copy'

        input :
        val assembly
        file gnomad_SNV_vcf
	file gnomad_SNV_index
	each chr

        output :
        file '*'

        script:
        if (${params.assembly}=="GRCh37" && $chr == "Y")
        """
        echo "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tAF\tAN\tAC\tnhomalt" > ${gnomad_SNV_vcf.simpleName}_frequency_table_Y.tsv
        echo "Y\t1\t.\tN\tN\t.\t.\t.\t.\t.\t." >> ${gnomad_SNV_vcf.simpleName}_frequency_table_Y.tsv
        """
        else
        """
        gatk --java-options "-Xmx4G" \
	VariantsToTable \
        -V ${gnomad_SNV_vcf} \
        -O ${gnomad_SNV_vcf.simpleName}_frequency_table_${chr}.tsv \
	--show-filtered \
	-L ${chr} \
	-F CHROM \
        -F POS \
        -F ID \
	-F REF \
        -F ALT \
        -F QUAL \
        -F FILTER \
        -F AF \
        -F AC \
        -F AN \
        -F nhomalt
	"""
}
