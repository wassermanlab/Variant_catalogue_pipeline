#R script
#Created by Solenne Correard in December 2021
#Owned by the Silent Genomes Project Activity 3 team
#Developped to build the IBVL, a background variant library

library(vcfR)
library(dplyr)
library(stringr)
library(tidyr)
options(scipen = 999)
###################
# Create the table as expected in the SQL
###################

#Read the arguments from the nextflow process calling the Rscript named 'MT_heteroplasmy_bin.nf'
args <- commandArgs(trailingOnly = TRUE)

#Define assembly from the args
assembly=(args[1])

####Organize different tables
##Read gnomAD MT table
gnomad_file=read.table(args[2], header=TRUE)
#Create the variant ID (chr-Pos_ref_alt)
ID_table_gnomad=gnomad_file[,c("chromosome", "position", "ref", "alt")]
ID_db_gnomad=paste(ID_table_gnomad$chromosome, ID_table_gnomad$position, ID_table_gnomad$ref, ID_table_gnomad$alt, sep="_")
gnomad_file=cbind(ID_db_gnomad, gnomad_file)

##MT annotation table (IBVL)
frequ_annot_file=read.vcfR(args[3])

#Severity table
severity_table=read.table((args[5]), fill=TRUE, header=TRUE)

table_frequ_db=data.frame()
table_annot_MT=data.frame()
variants_table=data.frame()
table_variant_transcript=data.frame()
table_variant_consequence=data.frame()
table_variant_annotation=data.frame()

vcf = frequ_annot_file@fix

chr = vcf[,"CHROM"]
pos = as.numeric(vcf[,"POS"])
variant = vcf[,"ID"]
ref = vcf[,"REF"]
alt = vcf[,"ALT"]

gt = frequ_annot_file@gt
AC_hom = gt[,"AC_hom"]
AC_het = gt[,"AC_het"]
AF_hom = gt[,"AF_hom"]
AF_het = gt[,"AF_het"]
AN = gt[,"AN"]
max_observed_heteroplasmy = gt[,"max_observed_heteroplasmy"]
heteroplasmy_histogram = gt[,"heteroplasmy_histogram"]

info = vcf[,"filters"]
filters = as.data.frame(sapply(strsplit(info, ";"),"[[",1))

##Generate table with info and annot
#VEP annotation
vep_annot = as.data.frame(sapply(strsplit(info, "CSQ"),"[[",2))

show("vep annot ok")

#Split annot as one varaint can have several annot
#Create large table
all_info = cbind(chr, pos, variant, ref, alt, filters, AC_hom, AC_het, AF_hom, AF_het, AN, max_observed_heteroplasmy, heteroplasmy_histogram, vep_annot)
colnames(all_info) = c("chr", "pos", "variant", "ref", "alt", "filters", "AC_hom", "AC_het", "AF_hom", "AF_het", "AN", "max_observed_heteroplasmy", "heteroplasmy_histogram", "vep_annot")

show("all_info ok")

#Creates one row per annotation
all_info_split = separate_rows(all_info, vep_annot, sep = ",")
vep_annot = all_info_split$vep_annot

#Split the VEP annotation (seperated by |)
vep_annot_split = do.call(rbind.data.frame, strsplit(vep_annot, "\\|"))
colnames(vep_annot_split) = c("Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type", "Feature", "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation", "DISTANCE", "STRAND", "FLAGS", "VARIANT_CLASS", "SYMBOL_SOURCE", "HGNC_ID", "TSL", "REFSEQ_MATCH", "SOURCE", "REFSEQ_OFFSET", "GIVEN_REF", "USED_REF", "BAM_EDIT", "SIFT", "PolyPhen", "HGVS_OFFSET", "CLIN_SIG", "SOMATIC", "PHENO", "VAR_SYNONYMS", "CADD_PHRED", "CADD_RAW", "SpliceAI_pred_DP_AG", "SpliceAI_pred_DP_AL", "SpliceAI_pred_DP_DG", "SpliceAI_pred_DP_DL", "SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL", "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL")
			      #, "SpliceAI_pred_SYMBOL"

show("vep_annot_split ok")

#Create large table
all_info_complete = cbind(all_info_split$chr, all_info_split$pos, all_info_split$variant, all_info_split$ref, all_info_split$alt, all_info_split$filters, all_info_split$AC_hom, all_info_split$AC_het, all_info_split$AF_hom, all_info_split$AF_het, all_info_split$AN, all_info_split$max_observed_heteroplasmy, all_info_split$heteroplasmy_histogram, vep_annot_split$Allele, vep_annot_split$Consequence, vep_annot_split$IMPACT, vep_annot_split$SYMBOL, vep_annot_split$Gene, vep_annot_split$Feature_type, vep_annot_split$Feature, vep_annot_split$BIOTYPE, vep_annot_split$EXON, vep_annot_split$INTRON, vep_annot_split$HGVSc, vep_annot_split$HGVSp, vep_annot_split$cDNA_position, vep_annot_split$CDS_position, vep_annot_split$Protein_position, vep_annot_split$Amino_acids, vep_annot_split$Codons, vep_annot_split$Existing_variation, vep_annot_split$DISTANCE, vep_annot_split$STRAND, vep_annot_split$FLAGS, vep_annot_split$VARIANT_CLASS, vep_annot_split$SYMBOL_SOURCE, vep_annot_split$HGNC_ID, vep_annot_split$TSL, vep_annot_split$REFSEQ_MATCH, vep_annot_split$SOURCE, vep_annot_split$REFSEQ_OFFSET, vep_annot_split$GIVEN_REF, vep_annot_split$USED_REF, vep_annot_split$BAM_EDIT, vep_annot_split$SIFT, vep_annot_split$PolyPhen, vep_annot_split$HGVS_OFFSET, vep_annot_split$CLIN_SIG, vep_annot_split$SOMATIC, vep_annot_split$PHENO, vep_annot_split$VAR_SYNONYMS, vep_annot_split$CADD_PHRED, vep_annot_split$CADD_RAW, vep_annot_split$SpliceAI_pred_DP_AG, vep_annot_split$SpliceAI_pred_DP_AL, vep_annot_split$SpliceAI_pred_DP_DG, vep_annot_split$SpliceAI_pred_DP_DL, vep_annot_split$SpliceAI_pred_DS_AG, vep_annot_split$SpliceAI_pred_DS_AL, vep_annot_split$SpliceAI_pred_DS_DG, vep_annot_split$SpliceAI_pred_DS_DL)

colnames(all_info_complete) = c("chr", "pos", "variant", "ref", "alt", "filters", "ac_hom", "ac_het", "af_hom", "af_het", "an", "max_hl", "heteroplasmy_histogram", "Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type", "Feature", "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation", "DISTANCE", "STRAND", "FLAGS", "VARIANT_CLASS", "SYMBOL_SOURCE", "HGNC_ID", "TSL", "REFSEQ_MATCH", "SOURCE", "REFSEQ_OFFSET", "GIVEN_REF", "USED_REF", "BAM_EDIT", "SIFT", "PolyPhen", "HGVS_OFFSET", "CLIN_SIG", "SOMATIC", "PHENO", "VAR_SYNONYMS", "CADD_PHRED", "CADD_RAW", "SpliceAI_pred_DP_AG", "SpliceAI_pred_DP_AL", "SpliceAI_pred_DP_DG", "SpliceAI_pred_DP_DL", "SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL", "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL")

All_info = as.data.frame(all_info_complete)

show("all info ok")


All_info$variant = paste(All_info$chr, All_info$pos, All_info$ref, All_info$alt, sep = "_")

All_info$pos = as.numeric(All_info$pos)

#In gnomAD, indels are identified as M-105-CGGAGCA-C
All_info$variant= with(All_info, ifelse((nchar(as.character(All_info$ref))>1 & nchar(as.character(All_info$alt))==1), paste("chrM", All_info$pos+1, substring(All_info$ref, 2), "-", sep="_"), 
					ifelse((nchar(as.character(All_info$ref))==1 & nchar(as.character(All_info$alt))>1), paste("chrM", All_info$pos+1, "-", substring(All_info$alt, 2), sep="_"), All_info$variant)))

All_info$hl_hist_temp = gsub("[\\[\\]]", "", regmatches(All_info$heteroplasmy_histogram, gregexpr("\\[.*?\\]", All_info$heteroplasmy_histogram))[[1]])[2]
All_info$hl_hist  = substring(All_info$hl_hist_temp, 2, nchar(All_info$hl_hist_temp)-1)

# Intergenic (Y/N) : From annotation file, if there is a cDNA_position --> Y, else N
All_info$intergenic = with(All_info, ifelse((max(as.character(All_info$cDNA_position)) > 0), "N", "Y"))

# UCSC URL : https://genome.ucsc.edu/cgi-bin/hgTracks?db=<assembly>&highlight=<assembly>.chrM%3A<pos>-<pos>&position=chrM%3A<pos-25>-<pos+25> / hg38 or hg19 db=hg38&highlight=hg38.chrM%3A8602-8602&position=chrM%3A8577-8627
All_info$ucsc_url=paste0("https://genome.ucsc.edu/cgi-bin/hgTracks?db=",assembly, "&highlight=", assembly, ".chrM%3A", All_info$pos, "-", All_info$pos, "&position=chrM%3A", All_info$pos-25, "-", All_info$pos+25)
  
# mitomap URL : https://mitomap.org/cgi-bin/search_allele?variant=<pos><ref>%3E<alt> // 8602T%3EC
All_info$mitomap_url=paste0("https://mitomap.org/cgi-bin/search_allele?variant=", All_info$pos, All_info$ref, "%3E", All_info$alt)
  
# gnomad_URL : https://gnomad.broadinstitute.org/variant/M-<pos>-<ref>-<alt>?dataset=gnomad_r3 // M-8602-T-C
#Only displayed if variant is in gnomad table
All_info$gnomad_url= with(All_info, ifelse((All_info$variant %in% gnomad_file$ID_db_gnomad), paste0("https://gnomad.broadinstitute.org/variant/M-", All_info$pos, "-", All_info$ref, "-", All_info$alt, "?dataset=gnomad_r3"), "NA"))

  	
# dbsnp_id : From annotation file, "Existing_variation" column
# dbsnp_URL : https://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=<rs_number> // rs1556423501
All_info$dbsnp_id= with(All_info, ifelse((grepl("rs", All_info$Existing_variation)), gsub(",.*$", "", All_info$Existing_variation), "NA"))
All_info$dbsnp_url= with(All_info, ifelse((grepl("rs", All_info$Existing_variation)), paste0("https://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=", All_info$dbsnp_id), "NA"))

#clinvar_VCV : From annotation file, "VCV" info
#If clinvar number is specified (If column contains ClinVar)
All_info$clinvar_vcv= with(All_info, ifelse((grepl("ClinVar::", All_info$VAR_SYNONYMS)), str_extract(All_info$VAR_SYNONYMS, "(?<=VCV)[0-9]*"), "NA"))
All_info$clinvar_url= with(All_info, ifelse((grepl("ClinVar::", All_info$VAR_SYNONYMS)),paste0("https://www.ncbi.nlm.nih.gov/clinvar/variation/", All_info$clinvar_vcv, "/"), "NA"))
	 
#transcript=All_info$Feature

#Replace empty cells by "NA"
All_info[All_info == ""] <- "NA"

#Create tables
# MT_IBVL_frequency
# Variant ID, AN, AC Hom, AC Het, AF Hom, AF het, hl hist, max hl
table_frequ_db = cbind(All_info$variant, All_info$an, All_info$ac_hom, All_info$ac_het, All_info$af_hom,  All_info$af_het, All_info$hl_hist, All_info$max_hl)
colnames(table_frequ_db) = c("variant", "an", "ac_hom", "ac_het", "af_hom", "af_het", "hl_hist", "max_hl")
table_frequ_db = unique(table_frequ_db)
table_frequ_db = as.data.frame(table_frequ_db)
#Remove the varaints with AN=0
table_frequ_db_AN0 = table_frequ_db[!grepl("0",table_frequ_db$an), ]
write.table(table_frequ_db_AN0, file="mt_ibvl_frequencies.tsv", quote=FALSE, row.names = FALSE, sep="\t")

show("frequ ok")

# MT_annotation
# variant_ID, pos, ref, alt, intergenic (Y/N), UCSC_URL, mitomap_URL, gnomad_URL, dbsnp_id, dbsnp_url, clinvar_url
table_annot_MT = cbind(All_info$variant, All_info$pos, All_info$ref, All_info$alt, All_info$ucsc_url, All_info$mitomap_url, All_info$gnomad_url, All_info$dbsnp_id, All_info$dbsnp_url, All_info$clinvar_url, All_info$clinvar_vcv)
colnames(table_annot_MT) = c("All_info$variant", "pos", "ref", "alt", "ucsc_url", "mitomap_url", "gnomad_url", "dbsnp_id", "dbsnp_url", "clinvar_url", "clinvar_vcv")
table_annot_MT = unique(table_annot_MT)
write.table(table_annot_MT, file="mts.tsv", quote=FALSE, row.names = FALSE, sep="\t")

show("annot ok")

#Variants Table
# Variant_id, variant_type (SNV/MT/SV)
variants_table=cbind(All_info$variant, "MT")
colnames(variants_table) = c("variant_id", "var_type")
variants_table = unique(variants_table)
write.table(variants_table, file="variants_MT.tsv", quote=FALSE, row.names = FALSE, sep="\t")

show("variants ok")

#Variants_transcript table
#transcript_id, variant_id, hgvsc
table_variant_transcript=cbind(All_info$Feature, All_info$variant, All_info$HGVSc)
colnames(table_variant_transcript) = c("transcript", "variant", "hgvsc")
table_variant_transcript = unique(table_variant_transcript)
table_variant_transcript = as.data.frame(table_variant_transcript)
table_variant_transcript_noNA = table_variant_transcript[!grepl("NA",table_variant_transcript$transcript), ]
write.table(table_variant_transcript_noNA, file="variants_transcripts_MT.tsv", quote=FALSE, row.names = FALSE, sep="\t")

show("var transcripts ok")

# Variants_consequences table
# variant_transcript_id, severity (i.e consequence coded in number)
#Only for non intergenic variants
#If there is several consequences on the same line (separrated by a coma), create one line per consequence
table_variant_consequence = cbind(All_info$Consequence, All_info$variant, All_info$Feature)
colnames(table_variant_consequence) = c("consequence", "variant", "transcript")
table_variant_consequence_split = as.data.frame(table_variant_consequence)
table_variant_consequence_split=separate_rows(table_variant_consequence_split, consequence, sep = "&")
#Replace each consequence by it's severity number
table_variant_severities <- table_variant_consequence_split
table_variant_severities$severity <- severity_table$severity_number[match(table_variant_severities$consequence, severity_table$consequence)]
#Keep only the wanted info
table_variant_severities=table_variant_severities[,c("severity", "variant", "transcript")]
table_variant_severities = unique(table_variant_severities)
table_variant_severities = as.data.frame(table_variant_severities)
table_variant_severities_noNA= table_variant_severities[!grepl("NA",table_variant_severities$transcript), ]
write.table(table_variant_severities_noNA, file="variants_consequences_MT.tsv", quote=FALSE, row.names = FALSE, sep="\t")

show("severity ok")

#Variants_annotation
#variants_transcript_id, hgvsp, polyphen, sift (for polyphen and sift, score and interpretation together)
table_variant_annotation=cbind(All_info$HGVSp, All_info$SIFT, All_info$PolyPhen, All_info$Feature, All_info$variant)
colnames(table_variant_annotation) = c("hgvsp", "sift", "polyphen", "transcript", "variant")
table_variant_annotation = unique(table_variant_annotation)
table_variant_annotation =as.data.frame(table_variant_annotation)
table_variant_annotation_noNA = table_variant_annotation[!grepl("NA",table_variant_annotation$hgvsp), ]
write.table(table_variant_annotation_noNA, file="variants_annotations_MT.tsv", quote=FALSE, row.names = FALSE, sep="\t")

show("annot ok")

# MT_gnomAD_frequency
# Variant_ID, AN, AC_hom, AC_het, AF_hom, AF_het, max_hl
#Extract only the variant that are present in the IBVL
#gnomad_intersect = gnomad_file[gnomad_file$ID_db  %in%  as.data.frame(AF_table_ID)$ID_db_IBVL, ]
gnomad_intersect = gnomad_file[gnomad_file$ID_db_gnomad  %in% All_info$variant, ]
#Keep only the wanted info
gnomad_intersect_mini=gnomad_intersect[,c("ID_db_gnomad", "AN", "AC_hom", "AC_het", "AF_hom", "AF_het", "max_observed_heteroplasmy")]
#rename the columns as expected in the SQL
colnames(gnomad_intersect_mini)=c("variant", "an", "ac_hom", "ac_het", "af_hom", "af_het", "max_hl")
gnomad_intersect_mini = unique(gnomad_intersect_mini)
write.table(gnomad_intersect_mini, file="mt_gnomad_frequencies.tsv", quote=FALSE, row.names = FALSE, sep="\t")

show("gnomad ok")

# Gene table
# Short name
gene_table=as.data.frame(All_info$SYMBOL)
colnames(gene_table)=c("short_name")
gene_table = unique(gene_table)
gene_table = as.data.frame(gene_table)
gene_table_noNA = gene_table[!grepl("NA",gene_table$short_name), ]
gene_table_noNA = as.data.frame(gene_table_noNA)
colnames(gene_table_noNA)=c("short_name")
write.table(gene_table_noNA, file="genes_MT.tsv", quote=FALSE, row.names = FALSE, sep="\t")

show("gene ok")

#Trancript Table
#(gene_id??), transcript_id, transcript type  (E for Ensembl, R for Refseq)
#The transcript need to be associated to it's gene
#Need to replace Ensembl by E and refseq by R to save space
#The Transcript Support Level (TSL) is a method to highlight the well-supported and poorly-supported transcript models for users.
transcript_table=as.data.frame(unique(cbind(All_info$Feature, All_info$SYMBOL, All_info$SOURCE, All_info$TSL)))
colnames(transcript_table)=c("transcript_id", "gene", "transcript_type", "tsl")
transcript_table=transcript_table %>% mutate(transcript_type = str_replace(transcript_type, "Ensembl", "E"))
transcript_table=transcript_table %>% mutate(transcript_type = str_replace(transcript_type, "Refseq", "R"))
transcript_table = unique(transcript_table)
transcript_table = as.data.frame(transcript_table)
transcript_table_noNA = transcript_table[!grepl("NA",transcript_table$transcript_id), ]
write.table(transcript_table_noNA, file="transcripts_MT.tsv", quote=FALSE, row.names = FALSE, sep="\t")

show("transcript")

