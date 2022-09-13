#R script
#Created by Solenne Correard in December 2021
#Owned by the Silent Genomes Project Activity 3 team
#Developped to build the IBVL, a background variant library

#library(ggplot2)
#library("ggplot2")
# install.packages('vcfR', repos='https://cloud.r-project.org/')
library(vcfR)
library(dplyr)
library(stringr)
library(tidyr)
###################
# Create the table as expected in the SQL
# MT_IBVL_frequency
# Variant ID, AN, AC Hom, AC Het, AF Hom, AF het, hl hist, max hl
# AN, in gnomAD : Total number of individual with high quality sequence at the position : How to do as it is not joint calling? Should NA be considered as absent but good qual or absent?
# in the IBVL, temporarily : AN = 2X number of individuals in the cohort (constant)
# Homoplsmaic AC : in gnomAD : Number of individuals with homoplasmic or near-homoplasmic variant (heteroplasmy level > 0.95)
# heteroplasmic AC : in gnomAD : Number of individuals with a variant at heteroplasmy level 0.1 - 0.95
# Max observed heteroplamy, in gnomAD : Maximum heteroplamy level observed across all individuals

###################
# Create the table as expected in the SQL
# MT_gnomAD_frequency
# Variant_ID, AN, AC_hom, AC_het, AF_hom, AF_het, max_hl
# Download the data from gnomAD
# https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1/vcf/genomes/gnomad.genomes.v3.1.sites.chrM.reduced_annotations.tsv

###################
# Create the table as expected in the SQL
# MT_annotation
# variant_ID, pos, red, alt, intergenic (Y/N), UCSC_URL, mitomap_URL, gnomad_URL, dbsnp_id, dbsnp_url, clinvar_url
# Intergenic (Y/N) : From annotation file, if there is a cDNA_position --> Y, else N
# UCSC URL : https://genome.ucsc.edu/cgi-bin/hgTracks?db=<assembly>&highlight=<assembly>.chrM%3A<pos>-<pos>&position=chrM%3A<pos-25>-<pos+25> / hg38 or hg19 db=hg38&highlight=hg38.chrM%3A8602-8602&position=chrM%3A8577-8627
# mitomap URL : https://mitomap.org/cgi-bin/search_allele?variant=<pos>%<ref>%<alt> // 8602T%3EC
# gnomad_URL : https://gnomad.broadinstitute.org/variant/M-<pos>-<ref>-<alt>?dataset=gnomad_r3 // M-8602-T-C
# dbsnp_id : From annotation file, "Existing_variation" column
# dbsnp_URL : https://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=<rs_number> // rs1556423501
# TO add to the SQL : clinvar_VCV : From annotation file, "VCV" info
# clinvar URL : https://www.ncbi.nlm.nih.gov/clinvar/variation/<VCV>/ // 692920



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

##MT vcf from the database (IBVL)
#Read the vcf file with MT calls
#MT_vcf=read.vcfR(args[3])
#Ectract the AF
#AF_table=extract.gt(MT_vcf,element = "AF")
#Temporarily : Replace NA by zero because no joint calling
#AF_table[is.na(AF_table)] <- 0
#Create the variant ID (chr-Pos_ref_alt)
#ID_table_IBVL=as.data.frame(MT_vcf@fix[,c("CHROM", "POS", "REF", "ALT")])
#ID_db_IBVL=paste(ID_table_IBVL$CHROM, ID_table_IBVL$POS, ID_table_IBVL$REF, ID_table_IBVL$ALT, sep="_")
#AF_table_ID=cbind(ID_db_IBVL, AF_table)

#Can take Hail table with alrady calculated info rather than vcf with individual genotype info
#remove # at the beginning of the first line
frequ_file_temp = sub("#", "", readLines(args[3]))
show(frequ_file_temp)
frequ_file=read.table(text=frequ_file_temp, header=TRUE)
#Create the variant ID (chr-Pos_ref_alt)
#ID_table_Hail=as.data.frame(Hail_MT_output[,c("chromosome", "position", "ref", "alt")])
#ID_db_Hail=paste(ID_table_Hail$chromosome, ID_table_Hail$position, ID_table_Hail$ref, ID_table_Hail$alt, sep="_")
#Hail_MT_output_ID=cbind(ID_db_Hail, Hail_MT_output)


##MT annotation table (IBVL)
### !!!! Need to remove the hashtag in front of the header line outside of R
#Needed to extract rsID and clinvar ID
MT_raw_annotaton_file=read.table((args[4]), fill=TRUE, header=TRUE)
#Change the Uploaded_variation format
MT_raw_annotaton_file$Uploaded_variation=str_replace(MT_raw_annotaton_file$Uploaded_variation, "/" , "_")

#Severity table
severity_table=read.table((args[6]), fill=TRUE, header=TRUE)

table_frequ_db=data.frame()
table_annot_MT=data.frame()
variants_table=data.frame()
table_variant_transcript=data.frame()
table_variant_consequence=data.frame()
table_variant_annotation=data.frame()

for (i in 1:nrow(frequ_file)) {
	#show(i)
  	#variant=AF_table_ID[c(i),1]
  	#AF_table_i = AF_table[c(i),]
  	#ID_table_IBVL_i=ID_table_IBVL[c(i),]

  	#Define variables specific to variant i
  	#pos=as.numeric(as.character(ID_table_IBVL_i$POS))
  	#ref=ID_table_IBVL_i$REF
  	#alt=ID_table_IBVL_i$ALT
  	
	frequ_file_i=frequ_file[i,]
	show("frequ_file_i")
	show(frequ_file_i)

	variant = frequ_file_i$ID

	show("variant")
	show(variant)

	chrom = frequ_file_i$CHROM
	pos = frequ_file_i$POS
	ref = frequ_file_i$REF
	alt = frequ_file_i$ALT

	variant = paste(chrom, pos, ref, alt, sep = "_")

  	#Several annotation per variant, need to subset the annotation table to extract the annotation corresponding to variant i
  	#For indel, need to remove the letter common between the ref and alt and add 1 to the pos to match the vep output
  	#In gnomAD, indels are identified as M-105-CGGAGCA-C
  	if (nchar(as.character(ref))>1 & nchar(as.character(alt))==1){
    		modified_variant_POS=pos+1
    		modified_variant_REF=substring(ref, 2)
    		modified_variant_ALT="-"
    		variant=paste("chrM", modified_variant_POS, modified_variant_REF, modified_variant_ALT, sep="_")
    	} else if (nchar(as.character(ref))==1 & nchar(as.character(alt))>1) {
    		modified_variant_POS=pos+1
    		modified_variant_REF="-"
    		modified_variant_ALT=substring(alt, 2)
    		variant=paste("chrM", modified_variant_POS, modified_variant_REF, modified_variant_ALT, sep="_")
    	} 
	

	MT_raw_annotaton_i= MT_raw_annotaton_file[(MT_raw_annotaton_file$Uploaded_variation == pos),]

	show("MT_raw_annotaton_i")
        show(MT_raw_annotaton_i)



	#gnomad filters : Given these large numbers of false positive calls for variants with V AF < 0.10, for the initial release we chose to exclude samples with mtCN < 50 and to report only variants with V AF 0.10 as we have greater confidence that such variants represent genuine heteroplasmies and not NUMT-derived false positives
  	# AN : number of cells with value
  	#an = length(AF_table_i[AF_table_i>=0])*2
  	#AC hom
  	#ac_hom = sum(AF_table_i > 0.95)
  	#AF hom
  	#af_hom = ac_hom/an
  	#AC het
  	#ac_het = sum(AF_table_i > 0.1 & AF_table_i < 0.95)
  	#AF het
  	#af_het = ac_het/an
  	
	an = frequ_file_i$AN
	ac_hom = frequ_file_i$AC_hom
	af_hom = frequ_file_i$AF_hom
	ac_het = frequ_file_i$AC_het
	af_het = frequ_file_i$AF_het
	max_hl = frequ_file_i$max_observed_heteroplasmy

	hl_hist_temp = gsub("[\\[\\]]", "", regmatches(frequ_file_i$heteroplasmy_histogram, gregexpr("\\[.*?\\]", frequ_file_i$heteroplasmy_histogram))[[1]])[2]
	hl_hist  = substring(hl_hist_temp, 2, nchar(hl_hist_temp)-1)

	if (ac_hom==0 & ac_het ==0) {
		next
	}else {
  		#hl hist
  		#Calculate the stat for each row (Each variant) with heteroplasmy level superior to 0
  		#freq_i = hist(as.numeric(AF_table_i)[as.numeric(AF_table_i) > 0.1 ], breaks=seq(0,1,by=0.1), plot=FALSE)
  		#Extract the count from each bin
  		#counts_i =as.data.frame(t(freq_i$counts))
  		#hl_hist=paste(counts_i$V1, counts_i$V2, counts_i$V3, counts_i$V4, counts_i$V5, counts_i$V6, counts_i$V7, counts_i$V8, counts_i$V9, counts_i$V10, sep=",")
  		#Max hl
  		#max_hl = max(AF_table_i)
  
  		# Intergenic (Y/N) : From annotation file, if there is a cDNA_position --> Y, else N
  		if (max(as.character(MT_raw_annotaton_i$cDNA_position)) > 0) { 
      			intergenic = "N"
  		} else {
      			intergenic = "Y"
  		}

  		# UCSC URL : https://genome.ucsc.edu/cgi-bin/hgTracks?db=<assembly>&highlight=<assembly>.chrM%3A<pos>-<pos>&position=chrM%3A<pos-25>-<pos+25> / hg38 or hg19 db=hg38&highlight=hg38.chrM%3A8602-8602&position=chrM%3A8577-8627
  		ucsc_url=paste0("https://genome.ucsc.edu/cgi-bin/hgTracks?db=",assembly, "&highlight=", assembly, ".chrM%3A", pos, "-", pos, "&position=chrM%3A", pos-25, "-", pos+25)
  
  		# mitomap URL : https://mitomap.org/cgi-bin/search_allele?variant=<pos><ref>%3E<alt> // 8602T%3EC
  		mitomap_url=paste0("https://mitomap.org/cgi-bin/search_allele?variant=", pos, ref, "%3E", alt)
  
  		# gnomad_URL : https://gnomad.broadinstitute.org/variant/M-<pos>-<ref>-<alt>?dataset=gnomad_r3 // M-8602-T-C
  		#Only displayed if variant is in gnomad table
  		if (variant %in% gnomad_file$ID_db_gnomad) {
    			gnomad_url=paste0("https://gnomad.broadinstitute.org/variant/M-", pos, "-", ref, "-", alt, "?dataset=gnomad_r3")
  		} else {
    			gnomad_url="NA"
  		}
  	
  		# dbsnp_id : From annotation file, "Existing_variation" column
  		if (grepl("rs", MT_raw_annotaton_i$Existing_variation)) {
    			dbsnp_id = gsub(",.*$", "", MT_raw_annotaton_i$Existing_variation)
    			# dbsnp_URL : https://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=<rs_number> // rs1556423501
    			dbsnp_url=paste0("https://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=", dbsnp_id)
  		} else {
    			dbsnp_id="NA"
    			dbsnp_url="NA"
  		}
  
  		#clinvar_VCV : From annotation file, "VCV" info
  		#If clinvar number is specified (If column contains ClinVar)
  		if (grepl("ClinVar::", MT_raw_annotaton_i$VAR_SYNONYMS)) {
    			clinvar_vcv = str_extract(MT_raw_annotaton_i$VAR_SYNONYMS, "(?<=VCV)[0-9]*")
    			# clinvar URL : https://www.ncbi.nlm.nih.gov/clinvar/variation/<VCV>/ // 692920
    			clinvar_url=paste0("https://www.ncbi.nlm.nih.gov/clinvar/variation/", clinvar_vcv, "/")
  		} else {
    			clinvar_vcv = "NA"
    			clinvar_url = "NA"
  		}
  
		#transcript_id
                transcript=MT_raw_annotaton_i$Feature
                #hgvsc
                hgvsc=MT_raw_annotaton_i$HGVSc
                #variant_transcript_id
                variant_transcript=paste0(variant, "_", transcript)
                
		#Consequence
                consequence_i=MT_raw_annotaton_i$Consequence	
		
                #hgvsp
                hgvsp=MT_raw_annotaton_i$HGVSp
                #polyphen_score
                polyphen=MT_raw_annotaton_i$PolyPhen
                #SIFT score
                sift=MT_raw_annotaton_i$SIFT
                
  		#Create tables
  		# MT_IBVL_frequency
  		# Variant ID, AN, AC Hom, AC Het, AF Hom, AF het, hl hist, max hl
  		temp_table_frequ_db_i = cbind(variant, an, ac_hom, ac_het, af_hom,  af_het, hl_hist, max_hl)
  		table_frequ_db=rbind.data.frame(table_frequ_db, temp_table_frequ_db_i)

  		# MT_annotation
  		# variant_ID, pos, ref, alt, intergenic (Y/N), UCSC_URL, mitomap_URL, gnomad_URL, dbsnp_id, dbsnp_url, clinvar_url
  		temp_table_annot_MT_i = cbind(variant, pos, ref, alt, ucsc_url, mitomap_url, gnomad_url, dbsnp_id, dbsnp_url, clinvar_url, clinvar_vcv)
    		table_annot_MT=unique(rbind.data.frame(table_annot_MT, temp_table_annot_MT_i))

		#Variants Table
        	# Variant_id, variant_type (SNV/MT/SV)
        	variants_table_i=cbind(variant, "MT")
		variants_table=rbind.data.frame(variants_table, variants_table_i)

		#Variants_transcript table
                #transcript_id, variant_id, hgvsc
                temp_table_variant_transcript_i=cbind(transcript, variant, hgvsc)
                table_variant_transcript=unique(rbind.data.frame(table_variant_transcript, temp_table_variant_transcript_i))
	
	        # Variants_consequences table
                # variant_transcript_id, severity (i.e consequence coded in number)
		#Only for non intergenic variants
                #If there is several consequences on the same line (separrated by a coma), create one line per consequence
                table_variant_consequence_i=cbind(consequence_i, variant, transcript)
                table_variant_consequence=unique(rbind.data.frame(table_variant_consequence, table_variant_consequence_i))
                table_variant_consequence_split=separate_rows(table_variant_consequence, consequence_i, sep = ",")
		
	        #Variants_annotation
                #variants_transcript_id, hgvsp, polyphen, sift (for polyphen and sift, score and interpretation together)
                table_variant_annotation_i=cbind(hgvsp, sift, polyphen, transcript, variant)
                table_variant_annotation=unique(rbind.data.frame(table_variant_annotation, table_variant_annotation_i))	
	}
}

#Write tables
write.table(table_frequ_db, file="mt_ibvl_frequencies.tsv", quote=FALSE, row.names = FALSE, sep="\t")

write.table(table_annot_MT, file="mts.tsv", quote=FALSE, row.names = FALSE, sep="\t")

colnames(variants_table)=c("variant_id", "var_type")
write.table(variants_table, file="variants_MT.tsv", quote=FALSE, row.names = FALSE, sep="\t")

table_variant_transcript=table_variant_transcript[!grepl("^-$", table_variant_transcript$transcript),]
write.table(table_variant_transcript, file="variants_transcripts_MT.tsv", quote=FALSE, row.names = FALSE, sep="\t")


table_variant_consequence_split=table_variant_consequence_split[!grepl("^intergenic_variant$", table_variant_consequence_split$consequence_i),]
#Replace each consequence by it's severity number
table_variant_severities <- table_variant_consequence_split
table_variant_severities$severity <- severity_table$severity_number[match(table_variant_severities$consequence_i, severity_table$consequence)]
#Keep only the wanted info
table_variant_severities=table_variant_severities[,c("severity", "variant", "transcript")]
write.table(table_variant_severities, file="variants_consequences_MT.tsv", quote=FALSE, row.names = FALSE, sep="\t")

table_variant_annotation=table_variant_annotation[!grepl("^-$", table_variant_annotation$hgvsp),]
write.table(table_variant_annotation, file="variants_annotations_MT.tsv", quote=FALSE, row.names = FALSE, sep="\t")


# MT_gnomAD_frequency
# Variant_ID, AN, AC_hom, AC_het, AF_hom, AF_het, max_hl
#Extract only the variant that are present in the IBVL
#gnomad_intersect = gnomad_file[gnomad_file$ID_db  %in%  as.data.frame(AF_table_ID)$ID_db_IBVL, ]
gnomad_intersect = gnomad_file[gnomad_file$ID_db_gnomad  %in% table_frequ_db$variant, ]
#Keep only the wanted info
gnomad_intersect_mini=gnomad_intersect[,c("ID_db_gnomad", "AN", "AC_hom", "AC_het", "AF_hom", "AF_het", "max_observed_heteroplasmy")]
#rename the columns as expected in the SQL
colnames(gnomad_intersect_mini)=c("variant", "an", "ac_hom", "ac_het", "af_hom", "af_het", "max_hl")
write.table(gnomad_intersect_mini, file="mt_gnomad_frequencies.tsv", quote=FALSE, row.names = FALSE, sep="\t")

# Gene table
# Short name
gene_table=as.data.frame(unique(MT_raw_annotaton_file[,c("SYMBOL")]))
colnames(gene_table)=c("short_name")
gene_table=as.data.frame(gene_table[!grepl("^-$", gene_table$short_name),])
colnames(gene_table)=c("short_name")
write.table(gene_table, file="genes_MT.tsv", quote=FALSE, row.names = FALSE, sep="\t")

#Trancript Table
#(gene_id??), transcript_id, transcript type  (E for Ensembl, R for Refseq)
#The transcript need to be associated to it's gene
#Need to replace Ensembl by E and refseq by R to save space
#The Transcript Support Level (TSL) is a method to highlight the well-supported and poorly-supported transcript models for users.
transcript_table=as.data.frame(unique(MT_raw_annotaton_file[,c("Feature", "SYMBOL", "SOURCE", "TSL")]))
colnames(transcript_table)=c("transcript_id", "gene", "transcript_type", "tsl")
transcript_table=transcript_table %>% mutate(transcript_type = str_replace(transcript_type, "Ensembl", "E"))
transcript_table=transcript_table %>% mutate(transcript_type = str_replace(transcript_type, "Refseq", "R"))
transcript_table=transcript_table[!grepl("^-$", transcript_table$transcript_id),]
write.table(transcript_table, file="transcripts_MT.tsv", quote=FALSE, row.names = FALSE, sep="\t")



