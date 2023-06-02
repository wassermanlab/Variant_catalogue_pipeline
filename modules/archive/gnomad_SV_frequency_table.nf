// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// Extract the useful information from the gnomAD table to make it smaller and load it in R
// Extract first 40th columns of the bed

process gnomad_SV_frequency_table {
	tag "${gnomad_SV_bed}"

        input :
        file gnomad_SV_bed
	each chr

        output :
        path '*_nohash.tsv'

        script :
        """
	gzip -cd $gnomad_SV_bed | cut -d\$"\t" -f 1-40 > temp_gnomad_SV_frequ_table.tsv
	awk '{ if (\$1 == "#chrom" || \$1 == $chr) { print } }' temp_gnomad_SV_frequ_table.tsv > gnomad_SV_frequ_table_${chr}.tsv
	sed 's/#chrom/chrom/g' gnomad_SV_frequ_table_${chr}.tsv > gnomad_SV_frequ_table_${chr}_nohash.tsv
	"""
}
