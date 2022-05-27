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
frequ_file=read.table(args[3], header=TRUE)
chromosome=frequ_file[1,c("chrom")]

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
severity_table=read.table((args[7]), fill=TRUE, header=TRUE)


#for (i in 1:nrow(SNV_vcf@fix)) {
##Loop to split the number of variants
slots_var=c(seq(0, nrow(frequ_file), by=5000), nrow(frequ_file))
  
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
		variant=frequ_file$rsid[i]
		frequ_variant = frequ_file[frequ_file$rsid==variant,]
		#To remove the multiallelic info ffrom the ID	
		#variant = gsub(";.*$", "", variant)
		#GT_table_i = GT_table[c(i),]
		SNV_annot_i = SNV_raw_annotaton_file[SNV_raw_annotaton_file$Uploaded_variation==variant,]
 
		show(variant)
		show(head(SNV_raw_annotaton_file$Uploaded_variation))

		#Define variables specific to variant i
		chr = frequ_file$chrom[i]
		pos = frequ_file$pos[i]
		ref = frequ_file$ref[i]
		alt = frequ_file$alt[i]
		#chr=SNV_vcf@fix[i,c("CHROM")]
  		#pos=as.numeric(SNV_vcf@fix[i,c("POS")])
  		#ref=SNV_vcf@fix[i,c("REF")]
  		#alt=SNV_vcf@fix[i,c("ALT")]

  		#Variant quality
  		#quality = SNV_vcf@fix[i,c("QUAL")]
  
		#Frequency from Hail file
		quality = frequ_variant$qual
		af_total = frequ_variant$af_tot
		ac_total = frequ_variant$ac_tot
		an_total = frequ_variant$an_tot
		hom_alt_total = frequ_variant$hom_alt_tot
		af_xx = frequ_variant$af_xx
		ac_xx = frequ_variant$ac_xx
		an_xx = frequ_variant$an_xx
		hom_alt_xx  = frequ_variant$hom_alt_xx
		af_xy = frequ_variant$af_xy
		ac_xy =frequ_variant$af_xy
		an_xy = frequ_variant$an_xy
		hom_alt_xy = frequ_variant$hom_alt_xy
	
  		# Variant ID, AF_tot, AF_XX, AF_XY, AC_tot, AC_XX, AC_XY, AN_tot, AN_XX, AN_XY, Hom_alt_tot, Hom_alt_XX, Hom_alt_XY
  		# AN_tot : number of 0/0, 0/1 and 1/1 genotypes (avoid counting the ./.)
  		#an_total = 2*(sum(GT_table_i == "0/0", na.rm=T) + sum(GT_table_i == "0/1", na.rm=T) + sum(GT_table_i == "1/1", na.rm=T)) 
  		#AC tot
  		#ac_total = sum(GT_table_i == "0/1", na.rm=T) + 2*sum(GT_table_i == "1/1", na.rm=T)
  		#AF tot = AC/AN
  		#af_total = ac_total/an_total
  		#Number of individus homozygotes for the alternative allele (1/1)
  		#hom_alt_total = sum(GT_table_i == "1/1", na.rm=T) 

  		#For XX individuals
  		#For now, make fake false with individuals and sex : 
  		#sex_table =  read.table("sample_sex.tsv", header=TRUE)
  		#Subset the GT_Table for XX individuals
  		#XX_Samples = sex_table[sex_table$Sex=="XX",1]
  		#XX_GT_table_i = GT_table_i[XX_Samples]
  		# AN_XX
  		#an_XX = 2*(sum(XX_GT_table_i == "0/0", na.rm=T) + sum(XX_GT_table_i == "0/1", na.rm=T) + sum(XX_GT_table_i == "1/1", na.rm=T)) 
  		#AC XX
  		#ac_XX = sum(XX_GT_table_i == "0/1", na.rm=T) + 2*sum(XX_GT_table_i == "1/1", na.rm=T)
  		#AF X = AC/AN
  		#af_XX = ac_XX/an_XX
  		#Number of individus homozygotes for the alternative allele (1/1)
  		#hom_alt_XX = sum(XX_GT_table_i == "1/1", na.rm=T) 
  
  		#For XY individuals
  		#Subset the GT_Table for XY individuals
  		#XY_Samples = sex_table[sex_table$Sex=="XY",1]
  		#XY_GT_table_i = GT_table_i[XY_Samples]
  		# AN_XY
  		#an_XY = 2*(sum(XY_GT_table_i == "0/0", na.rm=T) + sum(XY_GT_table_i == "0/1", na.rm=T) + sum(XY_GT_table_i == "1/1", na.rm=T)) 
  		#AC XY
  		#ac_XY = sum(XY_GT_table_i == "0/1", na.rm=T) + 2*sum(XY_GT_table_i == "1/1", na.rm=T)
  		#AF X = AC/AN
  		#af_XY = ac_XY/an_XY
  		#Number of individus homozygotes for the alternative allele (1/1)
  		#hom_alt_XY = sum(XY_GT_table_i == "1/1", na.rm=T) 


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
  
  		# CADD_score / interpr
  		# If you would like to apply a cutoff on deleteriousness, e.g. to identify potentially pathogenic variants, we would suggest to put a cutoff somewhere between 10 and 20. Maybe at 15, as this also happens to be the median value for all possible canonical splice site changes and non-synonymous variants in CADD v1.0. However, there is not a natural choice here -- it is always arbitrary. We therefore recommend integrating C-scores with other evidence and to rank your candidates for follow up rather than hard filtering.
		cadd_score=SNV_annot_i$CADD_PHRED
  		if (cadd_score <=15) {
			cadd_intr = "Tolerable"
		} else if (cadd_score > 15) {
			cadd_intr = "Damaging"
		}

		#splice_ai=SNV_annot_i$splice_ai
		splice_ai="0.999"

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
  		temp_table_frequ_db_i = cbind(variant, af_total, af_xx, af_xy, ac_total, ac_xx, ac_xy, an_total, an_xx, an_xy, hom_alt_total, hom_alt_xx, hom_alt_xy, quality)
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
colnames(gnomad_intersect_mini)=c("variant", "af_total", "ac_total", "an_total", "hom_alt_total")
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

