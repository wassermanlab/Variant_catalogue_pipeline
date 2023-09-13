#R script
#Created by Solenne Correard in December 2021
#Owned by the Silent Genomes Project Activity 3 team
#Developped to build the IBVL, a background variant library

library(vcfR)
library(dplyr)
library(stringr)
library(tidyr)
###################
# Create the table as expected in the SQL
###################

#Read the arguments from the nextflow process calling the Rscript named 'MT_heteroplasmy_bin.nf'
args <- commandArgs(trailingOnly = TRUE)
if (!require(glue)) {
  install.packages("glue")  # Install the package if not installed
  library(glue)             # Load the package
}
if (!require(data.table)) {
  install.packages("data.table")  # Install the package if not installed
  library(data.table)             # Load the package
}
#Define assembly from the args
assembly=(args[1])
chr = args[2]
####Organize different tables
##Read gnomAD SNV vcf
pat = glue("gnomad_frequency_table_(chr)?{chr}.tsv")
gfile = grep(perl = TRUE,pattern = pat,x = list.files(),value = T)
gnomad_file=fread(gfile)
if (length(gnomad_file)>11)gnomad_file[,12]=NULL
names(gnomad_file)
head(gnomad_file)
#Create the variant ID (chr-Pos_ref_alt)
ID_table_gnomad=gnomad_file[,c("CHROM", "POS", "REF", "ALT")]
ID_db_gnomad=paste(ID_table_gnomad$CHROM, ID_table_gnomad$POS, ID_table_gnomad$REF, ID_table_gnomad$ALT, sep="_")
gnomad_file=cbind(ID_db_gnomad, gnomad_file)

#File with frequencies calculated in Hail and annotation from VEP
frequ_annot_file=read.vcfR(args[3])
chromosome=frequ_annot_file@fix[1,1]

#Severity table
severity_table=read.table((args[4]), fill=TRUE, header=TRUE)

table_frequ_SNV=data.frame()
table_annot_SNV=data.frame()
table_variant_transcript=data.frame()
table_variant_consequence=data.frame()
table_variant_annotation=data.frame()
able_variant_annotation_i=data.frame()
	
vcf = frequ_annot_file@fix
chr = vcf[,"CHROM"]
pos = as.numeric(vcf[,"POS"])
variant = vcf[,"ID"]
ref = vcf[,"REF"]
alt = vcf[,"ALT"]
quality = vcf[,"QUAL"]

info = vcf[,"INFO"]

#Info from Hail
#Frequ_table
#Frequ_table
af_tot_xx_xy = (sapply(strsplit(info, "AF_tot_XX_XY="),"[[",2))
af_tot_xx_xy = (sapply(strsplit(af_tot_xx_xy, ";"),"[[",1))

ac_tot_xx_xy = (sapply(strsplit(info, "AC_tot_XX_XY="),"[[",2))
ac_tot_xx_xy = (sapply(strsplit(ac_tot_xx_xy, ";"),"[[",1))

an_tot_xx_xy = (sapply(strsplit(info, "AN_tot_XX_XY="),"[[",2))
an_tot_xx_xy = (sapply(strsplit(an_tot_xx_xy, ";"),"[[",1))

hom_tot_xx_xy = (sapply(strsplit(info, "hom_tot_XX_XY="),"[[",2))
hom_tot_xx_xy = (sapply(strsplit(hom_tot_xx_xy, ";"),"[[",1))

af_tot = as.data.frame(sapply(strsplit(af_tot_xx_xy, ","),"[[",1))
ac_tot = as.data.frame(sapply(strsplit(ac_tot_xx_xy, ","),"[[",1))
an_tot = as.data.frame(sapply(strsplit(an_tot_xx_xy, ","),"[[",1))
hom_tot = as.data.frame(sapply(strsplit(hom_tot_xx_xy, ","),"[[",1))

af_xx = as.data.frame(sapply(strsplit(af_tot_xx_xy, ","),"[[",2))
ac_xx = as.data.frame(sapply(strsplit(ac_tot_xx_xy, ","),"[[",2))
an_xx = as.data.frame(sapply(strsplit(an_tot_xx_xy, ","),"[[",2))
hom_xx = as.data.frame(sapply(strsplit(hom_tot_xx_xy, ","),"[[",2))

af_xy = as.data.frame(sapply(strsplit(af_tot_xx_xy, ","),"[[",3))
ac_xy = as.data.frame(sapply(strsplit(ac_tot_xx_xy, ","),"[[",3))
an_xy = as.data.frame(sapply(strsplit(an_tot_xx_xy, ","),"[[",3))
hom_xy = as.data.frame(sapply(strsplit(hom_tot_xx_xy, ","),"[[",3))

show("frequ ok")

##Generate table with info and annot
#VEP annotation
vep_annot = as.data.frame(sapply(strsplit(info, "CSQ"),"[[",2))

#Split annot as one varaint can have several annot
#Create large table
all_info = cbind(chr, pos, variant, ref, alt, quality, af_tot, ac_tot, an_tot, hom_tot, af_xx, ac_xx, an_xx, hom_xx, af_xy, ac_xy, an_xy, hom_xy, vep_annot)

show("all_info before colnames ok")
colnames(all_info) = c ("chr", "pos", "variant", "ref", "alt", "quality", "af_tot", "ac_tot", "an_tot", "hom_tot", "af_xx", "ac_xx", "an_xx", "hom_xx", "af_xy", "ac_xy", "an_xy", "hom_xy", "vep_annot")

show("all_info ok")

#Creates one row per annotation
all_info_split = separate_rows(all_info, vep_annot, sep = ",")
vep_annot = all_info_split$vep_annot

#Split the VEP annotation (seperated by |)
vep_annot_split = do.call(rbind.data.frame, strsplit(vep_annot, "\\|"))
colnames(vep_annot_split) = c("Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type", "Feature", "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation", "DISTANCE", "STRAND", "FLAGS", "VARIANT_CLASS", "SYMBOL_SOURCE", "HGNC_ID", "TSL", "REFSEQ_MATCH", "SOURCE", "REFSEQ_OFFSET", "GIVEN_REF", "USED_REF", "BAM_EDIT", "SIFT", "PolyPhen", "HGVS_OFFSET", "CLIN_SIG", "SOMATIC", "PHENO", "VAR_SYNONYMS", "CADD_PHRED", "CADD_RAW", "SpliceAI_pred_DP_AG", "SpliceAI_pred_DP_AL", "SpliceAI_pred_DP_DG", "SpliceAI_pred_DP_DL", "SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL", "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL")
        #, "SpliceAI_pred_SYMBOL")

show("vep_annot_split ok")

#Create large table
all_info_complete = cbind(all_info_split$chr, all_info_split$pos, all_info_split$variant, all_info_split$ref, all_info_split$alt, all_info_split$quality, all_info_split$af_tot, all_info_split$ac_tot, all_info_split$an_tot, all_info_split$hom_tot, all_info_split$af_xx, all_info_split$ac_xx, all_info_split$an_xx, all_info_split$hom_xx, all_info_split$af_xy, all_info_split$ac_xy, all_info_split$an_xy, all_info_split$hom_xy, vep_annot_split$Allele, vep_annot_split$Consequence, vep_annot_split$IMPACT, vep_annot_split$SYMBOL, vep_annot_split$Gene, vep_annot_split$Feature_type, vep_annot_split$Feature, vep_annot_split$BIOTYPE, vep_annot_split$EXON, vep_annot_split$INTRON, vep_annot_split$HGVSc, vep_annot_split$HGVSp, vep_annot_split$cDNA_position, vep_annot_split$CDS_position, vep_annot_split$Protein_position, vep_annot_split$Amino_acids, vep_annot_split$Codons, vep_annot_split$Existing_variation, vep_annot_split$DISTANCE, vep_annot_split$STRAND, vep_annot_split$FLAGS, vep_annot_split$VARIANT_CLASS, vep_annot_split$SYMBOL_SOURCE, vep_annot_split$HGNC_ID, vep_annot_split$TSL, vep_annot_split$REFSEQ_MATCH, vep_annot_split$SOURCE, vep_annot_split$REFSEQ_OFFSET, vep_annot_split$GIVEN_REF, vep_annot_split$USED_REF, vep_annot_split$BAM_EDIT, vep_annot_split$SIFT, vep_annot_split$PolyPhen, vep_annot_split$HGVS_OFFSET, vep_annot_split$CLIN_SIG, vep_annot_split$SOMATIC, vep_annot_split$PHENO, vep_annot_split$VAR_SYNONYMS, vep_annot_split$CADD_PHRED, vep_annot_split$CADD_RAW, vep_annot_split$SpliceAI_pred_DP_AG, vep_annot_split$SpliceAI_pred_DP_AL, vep_annot_split$SpliceAI_pred_DP_DG, vep_annot_split$SpliceAI_pred_DP_DL, vep_annot_split$SpliceAI_pred_DS_AG, vep_annot_split$SpliceAI_pred_DS_AL, vep_annot_split$SpliceAI_pred_DS_DG, vep_annot_split$SpliceAI_pred_DS_DL)
			  
colnames(all_info_complete) = c("chr", "pos", "variant", "ref", "alt", "quality", "af_tot", "ac_tot", "an_tot", "hom_tot", "af_xx", "ac_xx", "an_xx", "hom_xx", "af_xy", "ac_xy", "an_xy", "hom_xy", "Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type", "Feature", "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation", "DISTANCE", "STRAND", "FLAGS", "VARIANT_CLASS", "SYMBOL_SOURCE", "HGNC_ID", "TSL", "REFSEQ_MATCH", "SOURCE", "REFSEQ_OFFSET", "GIVEN_REF", "USED_REF", "BAM_EDIT", "SIFT", "PolyPhen", "HGVS_OFFSET", "CLIN_SIG", "SOMATIC", "PHENO", "VAR_SYNONYMS", "CADD_PHRED", "CADD_RAW", "SpliceAI_pred_DP_AG", "SpliceAI_pred_DP_AL", "SpliceAI_pred_DP_DG", "SpliceAI_pred_DP_DL", "SpliceAI_pred_DS_AG", "SpliceAI_pred_DS_AL", "SpliceAI_pred_DS_DG", "SpliceAI_pred_DS_DL")
        #, "SpliceAI_pred_SYMBOL")  

All_info = as.data.frame(all_info_complete)

show("all info ok")

# variant_ID, type, length, chr, pos, ref, alt, cadd_score, cadd_interpr, dbsnp_id, dbsnp_url, UCSC_url, ensembl_url, clinvar_url, gnomad_url
# type :SNV, ins or del
#Length
# CADD_score / interpr
# dbsnp_id : From annotation file, "Existing_variation" column
# dbsnp_URL : https://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=<rs_number> // rs1556423501
# UCSC URL : https://genome.ucsc.edu/cgi-bin/hgTracks?db=<assembly>&highlight=<assembly>.chrM%3A<pos>-<pos>&position=chrM%3A<pos-25>-<pos+25> / hg38 or hg19 db=hg38&highlight=hg38.chrM%3A8602-8602&position=chrM%3A8577-8627
# Ensembl_url : https://uswest.ensembl.org/Homo_sapiens/Location/View?r=<chr>%3A<pos-25>-<pos+25> // r=17%3A63992802-64038237
# clinvar URL : https://www.ncbi.nlm.nih.gov/clinvar/variation/<VCV>/ // 692920
# gnomad_URL : https://gnomad.broadinstitute.org/variant/M-<pos>-<ref>-<alt>?dataset=gnomad_r3 // M-8602-T-C

#Length and type
All_info$length <- with(All_info, ifelse(All_info$VARIANT_CLASS == "SNV", "1",
					ifelse(All_info$VARIANT_CLASS == "insertion" , length(All_info$alt)-1,
					       ifelse(All_info$VARIANT_CLASS == "deletion" , length(All_info$ref)-1,
						      ifelse(All_info$VARIANT_CLASS == "indel" , abs(length(All_info$alt)-length(All_info$ref)), "NA")))))
	
# CADD_score / interpr
# If you would like to apply a cutoff on deleteriousness, e.g. to identify potentially pathogenic variants, we would suggest to put a cutoff somewhere between 10 and 20. Maybe at 15, as this also happens to be the median value for all possible canonical splice site changes and non-synonymous variants in CADD v1.0. However, there is not a natural choice here -- it is always arbitrary. We therefore recommend integrating C-scores with other evidence and to rank your candidates for follow up rather than hard filtering.

All_info$cadd_intr <- with(All_info, ifelse(as.numeric(All_info$CADD_PHRED) <=15, "Tolerable", "Damaging"))

#SpliceAI output a serie of score, and the highest one should be reported
#Output VCF with SpliceAI predictions ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL included in the INFO column (see table below for details). Only SNVs and simple INDELs (REF or ALT is a single base) within genes are annotated. Variants in multiple genes have separate predictions for each gene.
#Example : CA|TTN|0.07|1.00|0.00|0.00|-7|-1|35|-29


All_info$SpliceAI_pred_DS_AG = as.numeric(All_info$SpliceAI_pred_DS_AG)
All_info$SpliceAI_pred_DS_AL = as.numeric(All_info$SpliceAI_pred_DS_AL) 
All_info$SpliceAI_pred_DS_DG = as.numeric(All_info$SpliceAI_pred_DS_DG)
All_info$SpliceAI_pred_DS_DL = as.numeric(All_info$SpliceAI_pred_DS_DL)

All_info$splice_ai <- with(All_info, ifelse(((All_info$SpliceAI_pred_DS_AG > All_info$SpliceAI_pred_DS_AL) & (All_info$SpliceAI_pred_DS_AG > All_info$SpliceAI_pred_DS_DG) & (All_info$SpliceAI_pred_DS_AG > All_info$SpliceAI_pred_DS_DL)), All_info$SpliceAI_pred_DS_AG,
					    ifelse(((All_info$SpliceAI_pred_DS_AL>All_info$SpliceAI_pred_DS_AG ) & (All_info$SpliceAI_pred_DS_AL > All_info$SpliceAI_pred_DS_DG) & (All_info$SpliceAI_pred_DS_AL > All_info$SpliceAI_pred_DS_DL)), All_info$SpliceAI_pred_DS_AL,
						   ifelse(((All_info$SpliceAI_pred_DS_DG > All_info$SpliceAI_pred_DS_AG) & (All_info$SpliceAI_pred_DS_DG > All_info$SpliceAI_pred_DS_AL) & (All_info$SpliceAI_pred_DS_DG > All_info$SpliceAI_pred_DS_DL)), All_info$SpliceAI_pred_DS_DG, 
							  ifelse(((All_info$SpliceAI_pred_DS_DL > All_info$SpliceAI_pred_DS_AG) & (All_info$SpliceAI_pred_DS_DL > All_info$SpliceAI_pred_DS_AL) & (All_info$SpliceAI_pred_DS_DL > All_info$SpliceAI_pred_DS_DG)), All_info$SpliceAI_pred_DS_DL, "NA")))))
					    
        
# dbsnp_id : From annotation file (SNV_annot_i), "Existing_variation" column
# dbsnp_URL : https://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=<rs_number> // rs1556423501
All_info$dbsnp_id<- with(All_info, ifelse(grepl("rs", All_info$Existing_variation), gsub(",.*$", "", All_info$Existing_variation), "NA"))
All_info$dbsnp_url <- with(All_info, ifelse(grepl("rs", All_info$Existing_variation), paste0("https://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=", All_info$dbsnp_id), "NA"))

# UCSC URL : https://genome.ucsc.edu/cgi-bin/hgTracks?db=<assembly>&highlight=<assembly>.chrM%3A<pos>-<pos>&position=chrM%3A<pos-25>-<pos+25> / hg38 or hg19 db=hg38&highlight=hg38.chrM%3A8602-8602&position=chrM%3A8577-8627
All_info$pos= as.numeric(All_info$pos)

All_info$ucsc_url = paste0("https://genome.ucsc.edu/cgi-bin/hgTracks?db=",assembly, "&highlight=", assembly, ".chr", All_info$chr, "%3A", All_info$pos, "-", All_info$pos, "&position=chr", All_info$chr, "%3A", All_info$pos-25, "-", All_info$pos+25)
	
# Ensembl_url : https://uswest.ensembl.org/Homo_sapiens/Location/View?r=<chr>%3A<pos-25>-<pos+25> // r=17%3A63992802-64038237
All_info$ensembl_url<- with(All_info, ifelse(assembly=="GRCh38", paste0("https://uswest.ensembl.org/Homo_sapiens/Location/View?r=", chr, "%3A", All_info$pos-25, "-", All_info$pos+25), paste0("https://grch37.ensembl.org/Homo_sapiens/Location/View?r=", All_info$chr, "%3A", All_info$pos-25, "-", All_info$pos+25)))
	    
#clinvar_VCV : From annotation file, "VCV" info
#If clinvar number is specified (If column contains ClinVar)
#clinvar URL : https://www.ncbi.nlm.nih.gov/clinvar/variation/<VCV>/ // 692920	
All_info$clinvar_vcv<- with(All_info, ifelse(grepl("ClinVar::VCV", All_info$VAR_SYNONYMS), str_extract(All_info$VAR_SYNONYMS, "(?<=VCV)[0-9]*"), "NA"))	
All_info$clinvar_url<- with(All_info, ifelse(grepl("ClinVar::VCV", All_info$VAR_SYNONYMS), paste0("https://www.ncbi.nlm.nih.gov/clinvar/variation/", All_info$clinvar_vcv, "/"), "NA"))
	
# gnomad_URL : https://gnomad.broadinstitute.org/variant/M-<pos>-<ref>-<alt>?dataset=gnomad_r3 // M-8602-T-C
#Only displayed if variant is in gnomad table
All_info$gnomad_url <- with(All_info, ifelse(All_info$variant %in% gnomad_file$ID_db_gnomad & assembly=="GRCh38", paste0("https://gnomad.broadinstitute.org/variant/", All_info$chr, "-", All_info$pos, "-", All_info$ref, "-", All_info$alt, "?dataset=gnomad_r3"),
				   ifelse(All_info$variant %in% gnomad_file$ID_db_gnomad & assembly=="GRCh37", paste0("https://gnomad.broadinstitute.org/variant/", All_info$chr, "-", All_info$pos, "-", All_info$ref, "-", All_info$alt, "?dataset=gnomad_r2_1"), "NA")))

#Replace empty cells by "NA"
All_info[All_info == ""] <- "NA"  
	
#All_info$variant_transcript=paste0(All_info$variant, "_", All_info$Feature)


### Create tables
# SNV_IBVL_frequency
# Variant ID, AF_tot, AF_XX, AF_XY, AC_tot, AC_XX, AC_XY, AN_tot, AN_XX, AN_XY, Hom_alt_tot, Hom_alt_XX, Hom_alt_XY, qual
table_frequ_SNV = cbind(All_info$variant, All_info$af_tot, All_info$af_xx, All_info$af_xy, All_info$ac_tot, All_info$ac_xx, All_info$ac_xy, All_info$an_tot, All_info$an_xx, All_info$an_xy, All_info$hom_tot, All_info$hom_xx, All_info$hom_xy, All_info$quality)
colnames(table_frequ_SNV) = c("variant", "af_tot", "af_xx", "af_xy", "ac_tot", "ac_xx", "ac_xy", "an_tot", "an_xx", "an_xy", "hom_tot", "hom_xx", "hom_xy", "quality")
table_frequ_SNV = unique(table_frequ_SNV)
write.table(table_frequ_SNV, file=paste0("genomic_ibvl_frequencies_", chromosome, ".tsv"), quote=FALSE, row.names = FALSE, sep="\t")

show("IBVL table ok")

# SNV_annotation
# variant_ID, type, length, chr, pos, ref, alt, cadd_score, cadd_interpr, dbsnp_id, dbsnp_url, UCSC_url, ensembl_url, clinvar_url, gnomad_url
table_annot_SNV= cbind(All_info$variant, All_info$VARIANT_CLASS, All_info$length, All_info$chr, All_info$pos, All_info$ref, All_info$alt, All_info$CADD_PHRED, All_info$cadd_intr, All_info$dbsnp_id, All_info$dbsnp_url, All_info$ucsc_url, All_info$ensembl_url, All_info$clinvar_url, All_info$gnomad_url, All_info$clinvar_vcv, All_info$splice_ai)
colnames(table_annot_SNV) = c("variant", "type", "length", "chr", "pos", "ref", "alt", "cadd_score", "cadd_intr", "dbsnp_id", "dbsnp_url", "ucsc_url", "ensembl_url", "clinvar_url", "gnomad_url", "clinvar_vcv", "splice_ai")
table_annot_SNV = unique(table_annot_SNV)
write.table(table_annot_SNV, file=paste0("snvs_", chromosome,".tsv"), quote=FALSE, row.names = FALSE, sep="\t")

show("SNV annot ok")

#Variants_transcript table
#transcript_id, variant, hgvsc
#Remove rows without a transcript (intergenic variants)
table_variant_transcript= cbind(All_info$Feature, All_info$variant, All_info$HGVSc)
colnames(table_variant_transcript) = c("transcript", "variant", "hgvsc")
table_variant_transcript = as.data.frame (table_variant_transcript)
table_variant_transcript_noNA = table_variant_transcript[!grepl("NA", table_variant_transcript$transcript), ]
table_variant_transcript_noNA = unique(table_variant_transcript_noNA)
write.table(table_variant_transcript_noNA, file=paste0("variants_transcripts_", chromosome,".tsv"), quote=FALSE, row.names = FALSE, sep="\t")

show("Variant transcript ok")

# Variants_consequences table
# variant_transcript_id, severity (i.e consequence coded in number)
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
table_variant_severities_noNA = table_variant_severities[!grepl("NA",table_variant_severities$transcript), ]
write.table(table_variant_severities_noNA, file=paste0("variants_consequences_", chromosome,".tsv"), quote=FALSE, row.names = FALSE, sep="\t")

show("varaint consequence ok")

#Variants_annotation
#variants_transcript_id, hgvsp, polyphen, sift (for polyphen and sift, score and interpretation together)
table_variant_annotation =cbind(All_info$HGVSp, All_info$SIFT, All_info$PolyPhen, All_info$Feature, All_info$variant)
colnames(table_variant_annotation) = c("hgvsp", "sift", "polyphen", "transcript", "variant")
table_variant_annotation = as.data.frame(table_variant_annotation)
table_variant_annotation_noNA = table_variant_annotation[!grepl("NA",table_variant_annotation$hgvsp), ]
table_variant_annotation_noNA = unique(table_variant_annotation_noNA)
write.table(table_variant_annotation_noNA, file=paste0("variants_annotations_", chromosome,".tsv"), quote=FALSE, row.names = FALSE, sep="\t")

show("variant annotation ok")

#Variants
variants_table=cbind(as.data.frame(unique(All_info$variant)), "SNV")
colnames(variants_table)=c("variant_id", "var_type")
variants_table = unique(variants_table)
write.table(variants_table, file=paste0("variants_", chromosome,".tsv"), quote=FALSE, row.names = FALSE, sep="\t")

show("variant ok")

# SNV_gnomAD_frequency
# Variant_ID, AF_tot, AC_tot, AN_tot, Hom_alt_tot
#Extract only the variant that are present in the IBVL
#gnomad_intersect = gnomad_file[gnomad_file$ID_db_gnomad  %in%  SNV_vcf@fix[,c("ID")], ]
gnomad_intersect = gnomad_file[gnomad_file$ID_db_gnomad  %in%  All_info$variant, ]
#Keep only the wanted info
gnomad_intersect_mini=gnomad_intersect[,c("ID_db_gnomad", "AF", "AC", "AN", "nhomalt")]
#rename the columns as expected in the SQL
colnames(gnomad_intersect_mini)=c("variant", "af_tot", "ac_tot", "an_tot", "hom_tot")
gnomad_intersect_mini = unique(gnomad_intersect_mini)
write.table(gnomad_intersect_mini, file=paste0("genomic_gnomad_frequencies_", chromosome, ".tsv"), quote=FALSE, row.names = FALSE, sep="\t")

show("gnomad frequ ok")

# Gene table
# Short name, (NCBI gene Id, may not be necessary anymore)
gene_table=as.data.frame(unique(All_info$SYMBOL))
colnames(gene_table)=c("short_name")
gene_table = as.data.frame(gene_table)
gene_table_noNA = gene_table[!grepl("NA", gene_table$short_name),]
gene_table_noNA = unique(gene_table_noNA)
gene_table_noNA = as.data.frame(gene_table_noNA)
colnames(gene_table_noNA)=c("short_name")
write.table(gene_table_noNA, file=paste0("genes_", chromosome, ".tsv"), quote=FALSE, row.names = FALSE, sep="\t")

show("gene table ok")

#Trancript Table
#(gene_id??), transcript_id, transcript type  (E for Ensembl, R for Refseq)
#The transcript need to be associated to it's gene
#Need to replace Ensembl by E and refseq by R to save space
transcript_table=as.data.frame(unique(cbind(All_info$Feature, All_info$SYMBOL, All_info$SOURCE, All_info$TSL)))
colnames(transcript_table)=c("transcript_id", "gene", "transcript_type", "tsl")
transcript_table=transcript_table %>% mutate(transcript_type = str_replace(transcript_type, "Ensembl", "E"))
transcript_table=transcript_table %>% mutate(transcript_type = str_replace(transcript_type, "Refseq", "R"))
transcript_table = unique(transcript_table)
transcript_table = as.data.frame(transcript_table)
transcript_table_noNA = transcript_table[!grepl("NA", transcript_table$transcript_id),]
write.table(transcript_table_noNA, file=paste0("transcripts_", chromosome, ".tsv"), quote=FALSE, row.names = FALSE, sep="\t")

show("transcript table ok")

