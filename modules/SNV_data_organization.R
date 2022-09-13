#R script
#Created by Solenne Correard in December 2021
#Owned by the Silent Genomes Project Activity 3 team
#Developped to build the IBVL, a background variant library

#library(ggplot2)
#library("ggplot2")
#install.packages('tidyr', repos='https://cloud.r-project.org/')
# install.packages('vcfR', repos='https://cloud.r-project.org/')
library(vcfR)
library(dplyr)
library(stringr)
library(tidyr)
###################
# Create the table as expected in the SQL
# SNV_IBVL_frequency
# Variant ID, AF_tot, AF_XX, AF_XY, AC_tot, AC_XX, AC_XY, AN_tot, AN_XX, AN_XY, Hom_alt_tot, Hom_alt_XX, Hom_alt_XY, qual
# How to define XX and XY
# In gnomAD : https://broadinstitute.github.io/gnomad_methods/_modules/gnomad/sample_qc/sex.html
# Can use bcftools --guess-ploidy plugin (based on varaint calls)
#https://samtools.github.io/bcftools/howtos/install.html
# using XYalign (based on bam)
# https://github.com/SexChrLab/XYalign


###################
# Create the table as expected in the SQL
# SNV_gnomAD_frequency
# Variant_ID, AF_tot, AC_tot, AN_tot, Hom_alt_tot
# /mnt/common/DATABASES/REFERENCES/GRCh37/GNOMAD/V2.1.1/gnomad.genomes.r2.1.1.sites.vcf.gz
# Create a subset with only the wanted info (outside of R) --> R will load faster a smaller file

###################
# Create the table as expected in the SQL
# SNV_annotation
# variant_ID, type, length, chr, pos, ref, alt, cadd_score, cadd_interpr, dbsnp_id, dbsnp_url, UCSC_url, ensembl_url, clinvar_url, gnomad_url
# type :SNV, ins or del
# CADD_score / interpr
# dbsnp_id : From annotation file, "Existing_variation" column
# dbsnp_URL : https://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=<rs_number> // rs1556423501
# UCSC URL : https://genome.ucsc.edu/cgi-bin/hgTracks?db=<assembly>&highlight=<assembly>.chrM%3A<pos>-<pos>&position=chrM%3A<pos-25>-<pos+25> / hg38 or hg19 db=hg38&highlight=hg38.chrM%3A8602-8602&position=chrM%3A8577-8627
# Ensembl_url : https://uswest.ensembl.org/Homo_sapiens/Location/View?r=<chr>%3A<pos-25>-<pos+25> // r=17%3A63992802-64038237
# clinvar URL : https://www.ncbi.nlm.nih.gov/clinvar/variation/<VCV>/ // 692920
# gnomad_URL : https://gnomad.broadinstitute.org/variant/M-<pos>-<ref>-<alt>?dataset=gnomad_r3 // M-8602-T-C


#Read the arguments from the nextflow process calling the Rscript named 'MT_heteroplasmy_bin.nf'
args <- commandArgs(trailingOnly = TRUE)

#Define assembly from the args
#assembly="GRCh37"
assembly=(args[1])

####Organize different tables
##Read gnomAD SNV vcf
#gnomad_file=read.table("work/31/3153f698d3ef271dc40359ebf12440/gnomad_frequency_table.tsv", header=TRUE)
gnomad_file=read.table(args[2], header=TRUE)
#Create the variant ID (chr-Pos_ref_alt)
ID_table_gnomad=gnomad_file[,c("CHROM", "POS", "REF", "ALT")]
ID_db_gnomad=paste(ID_table_gnomad$CHROM, ID_table_gnomad$POS, ID_table_gnomad$REF, ID_table_gnomad$ALT, sep="_")
gnomad_file=cbind(ID_db_gnomad, gnomad_file)

#File with frequencies calculated in Hail
#frequ_file=read.table('vcf_to_try_hail/subset_SNV_mt_var_filtered_tot_XX_XY_info.tsv', header=TRUE)
frequ_file=read.vcfR(args[3])
#chromosome=frequ_file[1,c("chrom")]

##MT vcf from the database (IBVL)
#Read the vcf file with MT calls
#SNV_vcf=read.vcfR("/mnt/scratch/SILENT/Act3/Processed/Individual/GRCh37/Batch_DryRun/Run_20211220/DeepVariant/DeepVariant_GLnexus_Run_20211220.vcf.gz")
#SNV_vcf=read.vcfR(args[3])
#Ectract the GT
#GT_table=extract.gt(SNV_vcf,element = "GT")
#colnames(GT_table)= sub("_.*","",colnames(GT_table))
#GT_table_ID=cbind(SNV_vcf@fix[,c("ID")], GT_table)
#chromosome=SNV_vcf@fix[1,c("CHROM")]

##SNV annotation table (IBVL)
### !!!! Need to remove the hashtag in front of the header line outside of R
#Needed to extract rsID and clinvar ID
SNV_raw_annotaton_file=read.table(args[4], fill=TRUE, header=TRUE)
#SNV_raw_annotaton_file=read.table("work/47/1c2e7dd036d2ea5fe303cfeeeb8e3f/DeepVariant_GLnexus_Run_20211220_annotation_table_nohash.tsv", fill=TRUE, header=TRUE)

#sex_table=read.table("sample_sex.tsv", header=TRUE)
#sex_table = read.table(args[5], header=TRUE)

#Severity table
severity_table=read.table((args[5]), fill=TRUE, header=TRUE)


#for (i in 1:nrow(SNV_vcf@fix)) {
##Loop to split the number of variants
slots_var=c(seq(0, nrow(SNV_raw_annotaton_file), by=5000), nrow(SNV_raw_annotaton_file))
  
for (j in 1:(length(slots_var)-1)){
    	#show(nrow(SNV_vcf@fix))
	#show(length(slots_var))
	#show("slot")
  	#show(j)
	min_i= slots_var[j]+1
	max_i= slots_var[j+1]

	table_frequ_SNV=data.frame()
        table_annot_SNV=data.frame()
	table_variant_transcript=data.frame()
	table_variant_consequence=data.frame()
        table_variant_annotation=data.frame()
	table_variant_annotation_i=data.frame()
	
	for (i in min_i: max_i){
		#show(i)
		SNV_annot_i = SNV_raw_annotaton_file[i,]
                #show("SNV_annot_i")
                #show(SNV_annot_i)

		variant = SNV_annot_i$Uploaded_variation

		frequ_file_i=frequ_file@fix[frequ_file@fix[,c("ID")]==variant,]
		#show("frequ_file_i")
                #show(frequ_file_i)


		#variant=frequ_file_i[c("ID")]
		#frequ_variant = frequ_file[frequ_file$rsid==variant,]
		#To remove the multiallelic info ffrom the ID	
		#variant = gsub(";.*$", "", variant)
		#GT_table_i = GT_table[c(i),]
		#SNV_annot_i = SNV_raw_annotaton_file[SNV_raw_annotaton_file$Uploaded_variation==variant,]
 

		#Define variables specific to variant i
		chr = frequ_file_i[c("CHROM")]
		pos = as.numeric(frequ_file_i[c("POS")])
		ref = frequ_file_i[c("REF")]
		alt = frequ_file_i[c("ALT")]

                # ????Intergenic ????
                #Type - VARIANT_CLASS
                type=SNV_annot_i$VARIANT_CLASS
		
		#length
                if (type=="SNV") {
                        length="1"
                } else if (type=="insertion") {
                        length = length(alt)-1
                } else if (type=="deletion") {
                        length = length(ref)-1
                }

		#If the varaint is longer than 49bp, then, it will be classified as a SV and should not be called by the SNV pipeline
		#Varaints longer than 49 should have been removed by the hail step already
		if (length > 49){
			next
		} else {

			#Define variables specific to variant i
			#Frequency from Hail file

			quality = frequ_file_i[c("QUAL")]
			info_frequ=frequ_file_i[c("INFO")]
			AC_i = unlist(strsplit(info_frequ, ";"))[1]
			AF_i = unlist(strsplit(info_frequ, ";"))[2]
			AN_i = unlist(strsplit(info_frequ, ";"))[3]
			hom_i = unlist(strsplit(info_frequ, ";"))[4]

			af_tot_temp = unlist(strsplit(AF_i, ","))[1]
			af_tot = unlist(strsplit(af_tot_temp, "="))[2]
			ac_tot_temp = unlist(strsplit(AC_i, ","))[1]
			ac_tot = unlist(strsplit(ac_tot_temp, "="))[2]
			an_tot_temp = unlist(strsplit(AN_i, ","))[1]
			an_tot = unlist(strsplit(an_tot_temp, "="))[2]
			hom_tot_temp = unlist(strsplit(hom_i, ","))[1]
			hom_tot = unlist(strsplit(hom_tot_temp, "="))[2]

			af_xx = unlist(strsplit(AF_i, ","))[2]
			ac_xx = unlist(strsplit(AC_i, ","))[2]
			an_xx = unlist(strsplit(AN_i, ","))[2]
			hom_xx = unlist(strsplit(hom_i, ","))[2]

			af_xy = unlist(strsplit(AF_i, ","))[3]
			ac_xy = unlist(strsplit(AC_i, ","))[3]
			an_xy = unlist(strsplit(AN_i, ","))[3]
			hom_xy = unlist(strsplit(hom_i, ","))[3]
  	
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
    
  			# CADD_score / interpr
  			# If you would like to apply a cutoff on deleteriousness, e.g. to identify potentially pathogenic variants, we would suggest to put a cutoff somewhere between 10 and 20. Maybe at 15, as this also happens to be the median value for all possible canonical splice site changes and non-synonymous variants in CADD v1.0. However, there is not a natural choice here -- it is always arbitrary. We therefore recommend integrating C-scores with other evidence and to rank your candidates for follow up rather than hard filtering.
			
			cadd_score=SNV_annot_i$CADD_PHRED
			if (cadd_score <=15) {
				cadd_intr = "Tolerable"
			} else if (cadd_score > 15) {
				cadd_intr = "Damaging"
			}

			#SpliceAI output a serie of score, and the highest one should be reported
			#Output VCF with SpliceAI predictions ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL included in the INFO column (see table below for details). Only SNVs and simple INDELs (REF or ALT is a single base) within genes are annotated. Variants in multiple genes have separate predictions for each gene.
			#Example : CA|TTN|0.07|1.00|0.00|0.00|-7|-1|35|-29

			if (grepl("|", SNV_annot_i$SpliceAI_pred)) {
				splice_ai = max(unlist(strsplit(SNV_annot_i$SpliceAI_pred, "\\|"))[3:6])
                        } else {
                                splice_ai = "NA"
                        }


  			# dbsnp_id : From annotation file (SNV_annot_i), "Existing_variation" column
  			if (grepl("rs", SNV_annot_i$Existing_variation)) {
    				dbsnp_id = gsub(",.*$", "", SNV_annot_i$Existing_variation)
    				# dbsnp_URL : https://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=<rs_number> // rs1556423501
    				dbsnp_url=paste0("https://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=", dbsnp_id)
  			} else {
    				dbsnp_id="NA"
    				dbsnp_url="NA"
  			}
  
  			# UCSC URL : https://genome.ucsc.edu/cgi-bin/hgTracks?db=<assembly>&highlight=<assembly>.chrM%3A<pos>-<pos>&position=chrM%3A<pos-25>-<pos+25> / hg38 or hg19 db=hg38&highlight=hg38.chrM%3A8602-8602&position=chrM%3A8577-8627
			ucsc_url=paste0("https://genome.ucsc.edu/cgi-bin/hgTracks?db=",assembly, "&highlight=", assembly, ".chr", chr, "%3A", pos, "-", pos, "&position=chr", chr, "%3A", pos-25, "-", pos+25)
  
  			# Ensembl_url : https://uswest.ensembl.org/Homo_sapiens/Location/View?r=<chr>%3A<pos-25>-<pos+25> // r=17%3A63992802-64038237
  			if (assembly=="GRCh38") {
    				ensembl_url=paste0("https://uswest.ensembl.org/Homo_sapiens/Location/View?r=", chr, "%3A", pos-25, "-", pos+25)
  			} else if (assembly=="GRCh37") {
    				ensembl_url=paste0("https://grch37.ensembl.org/Homo_sapiens/Location/View?r=", chr, "%3A", pos-25, "-", pos+25)
  			}  
    
  			#clinvar_VCV : From annotation file, "VCV" info
  			#If clinvar number is specified (If column contains ClinVar)
  			if (grepl("ClinVar::VCV", SNV_annot_i$VAR_SYNONYMS)) {
    				clinvar_vcv = str_extract(SNV_annot_i$VAR_SYNONYMS, "(?<=VCV)[0-9]*") 
    				#clinvar URL : https://www.ncbi.nlm.nih.gov/clinvar/variation/<VCV>/ // 692920
    				clinvar_url=paste0("https://www.ncbi.nlm.nih.gov/clinvar/variation/", clinvar_vcv, "/")
  			} else {
    				clinvar_vcv = "NA"
    				clinvar_url = "NA"
  			}
    
  			# gnomad_URL : https://gnomad.broadinstitute.org/variant/M-<pos>-<ref>-<alt>?dataset=gnomad_r3 // M-8602-T-C
  			#Only displayed if variant is in gnomad table
  			if (variant %in% gnomad_file$ID_db_gnomad) {
    				if (assembly=="GRCh38") {
      					#https://gnomad.broadinstitute.org/variant/1-55051215-G-GA?dataset=gnomad_r3
      					gnomad_url=paste0("https://gnomad.broadinstitute.org/variant/", chr, "-", pos, "-", ref, "-", alt, "?dataset=gnomad_r3")
    				} else if (assembly=="GRCh37") {
      					#https://gnomad.broadinstitute.org/variant/1-55516888-G-GA?dataset=gnomad_r2_1
      					gnomad_url=paste0("https://gnomad.broadinstitute.org/variant/", chr, "-", pos, "-", ref, "-", alt, "?dataset=gnomad_r2_1")
   				}
  			} else {
    				gnomad_url="NA"
  			}
  	
  			#transcript_id
  			transcript=SNV_annot_i$Feature
			#hgvsc
  			hgvsc=SNV_annot_i$HGVSc
  			#variant_transcript_id
			variant_transcript=paste0(variant, "_", transcript)
  			#Consequence
			consequence_i=SNV_annot_i$Consequence
			#hgvsp
			hgvsp=SNV_annot_i$HGVSp
			#polyphen_score
			polyphen=SNV_annot_i$PolyPhen
			#SIFT score
			sift=SNV_annot_i$SIFT

  			### Create tables
  			# SNV_IBVL_frequency
  			# Variant ID, AF_tot, AF_XX, AF_XY, AC_tot, AC_XX, AC_XY, AN_tot, AN_XX, AN_XY, Hom_alt_tot, Hom_alt_XX, Hom_alt_XY, qual
  			temp_table_frequ_db_i = cbind(variant, af_tot, af_xx, af_xy, ac_tot, ac_xx, ac_xy, an_tot, an_xx, an_xy, hom_tot, hom_xx, hom_xy, quality)
			table_frequ_SNV=unique(rbind.data.frame(table_frequ_SNV, temp_table_frequ_db_i))
	
	  		# SNV_annotation
	  		# variant_ID, type, length, chr, pos, ref, alt, cadd_score, cadd_interpr, dbsnp_id, dbsnp_url, UCSC_url, ensembl_url, clinvar_url, gnomad_url
	  		temp_table_annot_SNV_i = cbind(variant, type, length, chr, pos, ref, alt, cadd_score, cadd_intr, dbsnp_id, dbsnp_url, ucsc_url, ensembl_url, clinvar_url, gnomad_url, clinvar_vcv, splice_ai)
			table_annot_SNV=unique(rbind.data.frame(table_annot_SNV, temp_table_annot_SNV_i))
	
  			#Variants_transcript table
  			#transcript_id, variant, hgvsc
  			temp_table_variant_transcript_i=cbind(transcript, variant, hgvsc)
			table_variant_transcript=unique(rbind.data.frame(table_variant_transcript, temp_table_variant_transcript_i))

  			# Variants_consequences table
			# variant_transcript_id, severity (i.e consequence coded in number)
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
       
	### Write table slots
    	write.table(table_frequ_SNV, file=paste0("table_frequ_SNV_slot", j), quote=FALSE, row.names = FALSE)
    	write.table(table_annot_SNV, file=paste0("table_annot_SNV_slot", j), quote=FALSE, row.names = FALSE)
	table_variant_transcript_nodash=table_variant_transcript[!grepl("^-$", table_variant_transcript$transcript),]
	write.table(table_variant_transcript_nodash, file=paste0("table_variant_transcript_slot", j), quote=FALSE, row.names = FALSE)
	write.table(table_variant_consequence_split, file=paste0("table_variant_consequence_slot", j), quote=FALSE, row.names = FALSE)
        table_variant_annotation_nodash=table_variant_annotation[!grepl("^-$", table_variant_annotation$hgvsp),]
	write.table(table_variant_annotation_nodash, file=paste0("table_variant_annotation_slot", j), quote=FALSE, row.names = FALSE)
}

### Write tables 
# SNV_IBVL_frequency
list_frequ_tables_slots <- list.files(pattern = paste0("table_frequ_SNV_slot"))
tables_frequ_slots=lapply(list_frequ_tables_slots, read.table, header=TRUE)
combined_tables_frequ_slots=do.call(rbind, tables_frequ_slots)
write.table(combined_tables_frequ_slots, file=paste0("genomic_ibvl_frequencies_", chromosome, ".tsv"), quote=FALSE, row.names = FALSE, sep="\t")
file.remove(list_frequ_tables_slots)

# SNV_annotation
list_annot_tables_slots <- list.files(pattern = paste0("table_annot_SNV_slot"))
tables_annot_slots=lapply(list_annot_tables_slots, read.table, header=TRUE)
combined_tables_annot_slots=do.call(rbind, tables_annot_slots)
write.table(combined_tables_annot_slots, file=paste0("snvs_", chromosome,".tsv"), quote=FALSE, row.names = FALSE, sep="\t")
file.remove(list_annot_tables_slots)

#Variants_transcript table
list_variant_transcript_tables_slots <- list.files(pattern = paste0("table_variant_transcript_slot"))
tables_variant_transcript_slots=lapply(list_variant_transcript_tables_slots, read.table, header=TRUE)
combined_tables_variant_transcript_slots=do.call(rbind, tables_variant_transcript_slots)
combined_tables_variant_transcript_slots=combined_tables_variant_transcript_slots[!grepl("^-$", combined_tables_variant_transcript_slots),]
write.table(combined_tables_variant_transcript_slots, file=paste0("variants_transcripts_", chromosome,".tsv"), quote=FALSE, row.names = FALSE, sep="\t")
file.remove(list_variant_transcript_tables_slots)

# Variants_consequences table
list_variant_consequence_tables_slots <- list.files(pattern = paste0("table_variant_consequence_slot"))
tables_variant_consequence_slots=lapply(list_variant_consequence_tables_slots, read.table, header=TRUE)
combined_tables_variant_consequence_slots=do.call(rbind, tables_variant_consequence_slots)
combined_tables_variant_consequence_slots=combined_tables_variant_consequence_slots[!grepl("^intergenic_variant$", combined_tables_variant_consequence_slots$consequence_i),]
#Replace each consequence by it's severity number
table_variant_severities <- combined_tables_variant_consequence_slots
table_variant_severities$severity <- severity_table$severity_number[match(table_variant_severities$consequence_i, severity_table$consequence)]
#Keep only the wanted info
table_variant_severities=table_variant_severities[,c("severity", "variant", "transcript")]
write.table(table_variant_severities, file=paste0("variants_consequences_", chromosome,".tsv"), quote=FALSE, row.names = FALSE, sep="\t")
file.remove(list_variant_consequence_tables_slots)

#Variants_annotation
list_variant_annotation_tables_slots <- list.files(pattern = paste0("table_variant_annotation_slot"))
tables_variant_annotation_slots=lapply(list_variant_annotation_tables_slots, read.table, header=TRUE)
combined_tables_variant_annotation_slots=do.call(rbind, tables_variant_annotation_slots)
write.table(combined_tables_variant_annotation_slots, file=paste0("variants_annotations_", chromosome,".tsv"), quote=FALSE, row.names = FALSE, sep="\t")
file.remove(list_variant_annotation_tables_slots)

#Variants
variants_table=cbind(as.data.frame(unique(SNV_raw_annotaton_file[,c("Uploaded_variation")])), "SNV")
colnames(variants_table)=c("variant_id", "var_type")
write.table(variants_table, file=paste0("variants_", chromosome,".tsv"), quote=FALSE, row.names = FALSE, sep="\t")

# SNV_gnomAD_frequency
# Variant_ID, AF_tot, AC_tot, AN_tot, Hom_alt_tot
#Extract only the variant that are present in the IBVL
#gnomad_intersect = gnomad_file[gnomad_file$ID_db_gnomad  %in%  SNV_vcf@fix[,c("ID")], ]
gnomad_intersect = gnomad_file[gnomad_file$ID_db_gnomad  %in%  frequ_file[,c("rsid")], ]
#Keep only the wanted info
gnomad_intersect_mini=gnomad_intersect[,c("ID_db_gnomad", "AF", "AC", "AN", "nhomalt")]
#rename the columns as expected in the SQL
colnames(gnomad_intersect_mini)=c("variant", "af_tot", "ac_tot", "an_tot", "hom_tot")
write.table(gnomad_intersect_mini, file=paste0("genomic_gnomad_frequencies_", chromosome, ".tsv"), quote=FALSE, row.names = FALSE, sep="\t")

# Gene table
# Short name, (NCBI gene Id, may not be necessary anymore)
gene_table=as.data.frame(unique(SNV_raw_annotaton_file[,c("SYMBOL")]))
colnames(gene_table)=c("short_name")
gene_table=as.data.frame(gene_table[!grepl("^-$", gene_table$short_name),])
colnames(gene_table)=c("short_name")
write.table(gene_table, file=paste0("genes_", chromosome, ".tsv"), quote=FALSE, row.names = FALSE, sep="\t")

#Trancript Table
#(gene_id??), transcript_id, transcript type  (E for Ensembl, R for Refseq)
#The transcript need to be associated to it's gene
#Need to replace Ensembl by E and refseq by R to save space
transcript_table=as.data.frame(unique(SNV_raw_annotaton_file[,c("Feature", "SYMBOL", "SOURCE", "TSL")]))
colnames(transcript_table)=c("transcript_id", "gene", "transcript_type", "tsl")
transcript_table=transcript_table %>% mutate(transcript_type = str_replace(transcript_type, "Ensembl", "E"))
transcript_table=transcript_table %>% mutate(transcript_type = str_replace(transcript_type, "Refseq", "R"))
write.table(transcript_table, file=paste0("transcripts_", chromosome, ".tsv"), quote=FALSE, row.names = FALSE, sep="\t")

