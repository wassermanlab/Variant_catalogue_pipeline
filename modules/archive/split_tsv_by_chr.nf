// Nextflow process
// Created by Solenne Correard in December 2021
// Owned by the Silent Genomes Project Activity 3 team
// Developped to build the IBVL, a background variant library

// Overview of the process goal and characteristics :
// SNV calling, split the vcf by chr
//	split the multiallelic variants as one per line

process split_tsv_by_chr {

	input :
	file SNV_frequ_tsv
	val assembly
	val batch
	val run

	output :
	path '*.tsv', emit : vcf_onechr	

	script :
	"""
	#Replace ":" by a tab
	#sed -i 's/:/\t/' $SNV_frequ_tsv > SNV_frequ_chr.tsv

	#Split by first column
	awk 'FNR==1{hdr=\$0;next} {if (!seen[\$3]++) print hdr>\$3"_frequ.tsv"; print>\$3"_frequ.tsv"}' $SNV_frequ_tsv
	"""
}
