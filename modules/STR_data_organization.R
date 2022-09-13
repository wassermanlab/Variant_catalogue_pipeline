#R script
#Created by Solenne Correard in December 2021
#Owned by the Silent Genomes Project Activity 3 team
#Developped to build the IBVL, a background variant library

#library(ggplot2)
#library("ggplot2")
# install.packages('vcfR', repos='https://cloud.r-project.org/')
# install.packages('jsonlite', repos='https://cloud.r-project.org/')
#install.packages("data.table", repos="https://cloud.r-project.org/")
#library("data.table")
library(jsonlite)
library(vcfR)
library(dplyr)
library(stringr)
library(tidyr)
###################
# Create the table as expected in the SQL

#Read the arguments from the nextflow process calling the Rscript named 'MT_heteroplasmy_bin.nf'
args <- commandArgs(trailingOnly = TRUE)

#Define assembly from the args
#assembly="GRCh37"
assembly=(args[1])
var_type=(args[5])


#Read the vcf file with from the IBVL
STR_vcf=read.vcfR(args[2])
#Ectract the allele size distribution (REPCN = Repeat)
CN_table=extract.gt(STR_vcf,element = "REPCN")
CN_table_ID=cbind(STR_vcf@fix[,c("ID")], CN_table)

##For STR, the information is read from the JSON variant catalogue file
STR_json=fromJSON(args[3])
STR_json=separate_rows(STR_json, ReferenceRegion, sep = ",")

table_sv_consequence=data.frame()
table_str=data.frame()
variants_table=data.frame()
gene_table=data.frame()

for (i in 1: nrow(STR_vcf)){
	#show(i)
	variant=STR_vcf@fix[i,c("ID")]
	variant = gsub(";.*$", "", variant)
	
	#Define variables specific to variant i
	chr=STR_vcf@fix[i,c("CHROM")]
 	start=as.numeric(STR_vcf@fix[i,c("POS")])
  	ref=STR_vcf@fix[i,c("REF")]
  	alt=STR_vcf@fix[i,c("ALT")]

	end=extract.info(STR_vcf, element = "END")[i]
	#Number of repeat units spanned by the repeat in the reference
	n_STR_ref=extract.info(STR_vcf, element = "REF")[i]
	#Repeat unit in the reference orientation
	RU=extract.info(STR_vcf, element = "RU")[i]
	# Repeat id from the repeat-specification file
	REPID=extract.info(STR_vcf, element = "REPID")[i]

	#JSON info for this variant
	variant_reference_region=paste0(chr, ":", start, "-", end)
	json_info_i=STR_json[STR_json$ReferenceRegion==variant_reference_region,]

	gene=json_info_i$LocusId
#	GeneRegion=json_info_i$GeneRegion

	#Link to gnomAD STR page
	#https://gnomad.broadinstitute.org/short-tandem-repeat/DMD?dataset=gnomad_r3
	gnomad_url=paste0("https://gnomad.broadinstitute.org/short-tandem-repeat/", gene, "?dataset=gnomad_r3")

	#Allele size distrribution
	#For the graph, in the format : (allele size 1: n1, allele size 2 : n2, etc)
        CN_table_i = CN_table[c(i),]
        ref_geno=paste0(n_STR_ref, "/", n_STR_ref)

	#Create data frame with the length
	CN_table_i = as.data.frame(CN_table_i)

	#Genotype table : For each individual, allele 1 and allele 2
	table_genotype = cbind(CN_table_i, CN_table_i)
	colnames(table_genotype) = c("allele1", "allele2")
	table_genotype$allele1 <- table_genotype$allele1 %>% replace_na(ref_geno)
	table_genotype$allele2 <- table_genotype$allele2 %>% replace_na(ref_geno)

	table_genotype$allele1 = gsub("/.*", "", table_genotype$allele1)
	table_genotype$allele2 = gsub(".*/", "", table_genotype$allele2)
	
	#Allele frequency distribution
	data_frame_A1=as.data.frame(table_genotype$allele1)
	colnames(data_frame_A1) = c("Allele")
	data_frame_A2=as.data.frame(table_genotype$allele2)
	colnames(data_frame_A2) = c("Allele")
	
	Allele_table= rbind(data_frame_A1, data_frame_A2)
	Allele_count = as.data.frame(table(Allele_table))

	min_n_repeat = min(as.numeric(Allele_table[,1]), na.rm=T)
	max_n_repeat = max(as.numeric(Allele_table[,1]), na.rm=T)

	allele_size_distrib=","
	for (j in 1:nrow(Allele_count)){
		allele_size_distrib_j=paste0(Allele_count[j,1], ":", Allele_count[j,2], ",")
		allele_size_distrib=paste0(allele_size_distrib, allele_size_distrib_j)
	}

	#Remove first and last comma
	allele_size_distrib=str_sub(allele_size_distrib,2,-2)


  	### Create tables
  	# str
	#variant, repeat unit, min_n_repeat, max_n_repeat, allele distribution
	temp_str_i=cbind(variant, RU, min_n_repeat, max_n_repeat, allele_size_distrib)
	table_str = rbind.data.frame(table_str, temp_str_i)


	#Variants Table
	# Variant_id, variant_type (SNV/MT/SV)
	variants_table_i=cbind(variant, "SV")
	variants_table=unique(rbind.data.frame(variants_table, variants_table_i))       

	#sv_consequences
	show("gene")
	show(gene)
	show(length(gene))
	if (length(gene) > 0) {
		sv_consequences_i=cbind(gene, variant, "NA")
	}
	table_sv_consequence=unique(rbind.data.frame(table_sv_consequence, sv_consequences_i[1,]))

	#gene table
	gene_table_i=gene
	gene_table=unique(rbind.data.frame(gene_table, gene_table_i))
}


### Write tables 
# str
colnames(table_str)=c("variant", "repeat_unit", "min_n_repeat", "max_n_repeat", "allele_distrib")
write.table(table_str, file="str.tsv", quote=FALSE, row.names = FALSE)

#Variants
colnames(variants_table) = c("variant_id", "var_type")
write.table(variants_table, file="variants_str.tsv", quote=FALSE, row.names = FALSE)

#sv_consequences
colnames(table_sv_consequence)=c("gene", "variant", "consequence")
write.table(table_sv_consequence, file="sv_consequences_str.tsv", quote=FALSE, row.names = FALSE)

#gene table
colnames(gene_table)=c("short_name")
write.table(gene_table, file="genes_str.tsv", quote=FALSE, row.names = FALSE)
