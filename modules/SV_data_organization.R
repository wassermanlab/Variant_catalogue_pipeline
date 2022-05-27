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
# SV_IBVL_frequency
# Variant ID, AF_tot, AF_XX, AF_XY, AC_tot, AC_XX, AC_XY, AN_tot, AN_XX, AN_XY, Hom_alt_tot, Hom_alt_XX, Hom_alt_XY, qual
# How to define XX and XY
# In gnomAD : https://broadinstitute.github.io/gnomad_methods/_modules/gnomad/sample_qc/sex.html
# Can use bcftools --guess-ploidy plugin (based on varaint calls)
#https://samtools.github.io/bcftools/howtos/install.html
# using XYalign (based on bam)
# https://github.com/SexChrLab/XYalign


###################
# Create the table as expected in the SQL
# SV_gnomAD_frequency
#Need to decide if we display SV frequencies from gnomAD as the SV may not be exactly the same and could lead to confusion
# Variant_ID, AF_tot, AC_tot, AN_tot, Hom_alt_tot
# Create a subset with only the wanted info (outside of R) --> R will load faster a smaller file

###################
# Create the table as expected in the SQL
# SV_annotation (svs)
# variant_ID, chr1, chr1_pos1 (start), chr1_po2 (end), type, length, algorithm, ucsc_url, gnomad_id, gnomad_url


#Read the arguments from the nextflow process calling the Rscript named 'MT_heteroplasmy_bin.nf'
args <- commandArgs(trailingOnly = TRUE)

#Define assembly from the args
#assembly="GRCh37"
assembly=(args[1])
var_type=(args[6])

####Organize different tables
##Read gnomAD SNV vcf
#gnomad_file=read.table("work/31/3153f698d3ef271dc40359ebf12440/gnomad_frequency_table.tsv", header=TRUE)
#gnomad_file=read.table(args[2], header=TRUE)
#Create the variant ID (chr-Pos_ref_alt)
#ID_table_gnomad=gnomad_file[,c("chrom", "start", "end", "svtype")]
#ID_db_gnomad=paste(ID_table_gnomad$chrom, ID_table_gnomad$start, ID_table_gnomad$end, ID_table_gnomad$svtype, sep="_")
#gnomad_file=cbind(ID_db_gnomad, gnomad_file)

##MT vcf from the database (IBVL)
#Read the vcf file with MT calls
#SNV_vcf=read.vcfR("/mnt/scratch/SILENT/Act3/Processed/Individual/GRCh37/Batch_DryRun/Run_20211220/DeepVariant/DeepVariant_GLnexus_Run_20211220.vcf.gz")
SV_vcf=read.vcfR(args[2])
#Ectract the GT
GT_table=extract.gt(SV_vcf,element = "GT")
colnames(GT_table)= sub("_.*","",colnames(GT_table))
GT_table_ID=cbind(SV_vcf@fix[,c("ID")], GT_table)
chromosome=SV_vcf@fix[1,c("CHROM")]

##SNV annotation table (IBVL)
### !!!! Need to remove the hashtag in front of the header line outside of R
#Needed to extract rsID and clinvar ID
SV_raw_annotaton_file=read.table(args[3], fill=TRUE, header=TRUE)
#SNV_raw_annotaton_file=read.table("work/47/1c2e7dd036d2ea5fe303cfeeeb8e3f/DeepVariant_GLnexus_Run_20211220_annotation_table_nohash.tsv", fill=TRUE, header=TRUE)

#sex_table=read.table("sample_sex.tsv", header=TRUE)
sex_table = read.table(args[4], header=TRUE)

#for (i in 1:nrow(SNV_vcf@fix)) {
##Loop to split the number of variants
slots_var=c(seq(0, nrow(SV_vcf@fix), by=5000), nrow(SV_vcf@fix))
  
for (j in 1:(length(slots_var)-1)){
    	#show(nrow(SNV_vcf@fix))
	#show(length(slots_var))
	#show("slot")
  	#show(j)
	min_i= slots_var[j]+1
	max_i= slots_var[j+1]

	table_ctx=data.frame()
	table_frequ_SV=data.frame()
        table_annot_SV=data.frame()
	table_sv_consequence=data.frame()

	for (i in min_i: max_i){
		#show(i)
		variant=SV_vcf@fix[i,c("ID")]
		#To remove the multiallelic info ffrom the ID	
		#variant = gsub(";.*$", "", variant)
		GT_table_i = GT_table[c(i),]
		SV_annot_i = SV_raw_annotaton_file[SV_raw_annotaton_file$Uploaded_variation==variant,]
  
		#Define variables specific to variant i
		chr=SV_vcf@fix[i,c("CHROM")]
  		pos=as.numeric(SV_vcf@fix[i,c("POS")])
  		ref=SV_vcf@fix[i,c("REF")]
  		alt=SV_vcf@fix[i,c("ALT")]

  		#Variant quality
  		quality = SV_vcf@fix[i,c("QUAL")]
 
                #Type - VARIANT_CLASS (from annotation)
                type=SV_annot_i$VARIANT_CLASS
                #Type from vcf
                type_vcf=extract.info(SV_vcf, element = "SVTYPE")[i]

                #length
                #Info available in the INFO part of the vcf
		length = extract.info(SV_vcf, element = "SVLEN")[i]

		#Other info specific to SV
		AVG_LEN  = as.numeric(extract.info(SV_vcf, element = "AVG_LEN")[i])
		AVG_START = as.numeric(extract.info(SV_vcf, element = "AVG_START")[i])
		AVG_END = as.numeric(extract.info(SV_vcf, element = "AVG_END")[i])
		#This method indicate Jasmine, which is not the case
		SVMETHOD = extract.info(SV_vcf, element = "SVMETHOD")[i]
		#This is to extract the actual method
		IDLIST = extract.info(SV_vcf, element = "IDLIST")[i]
		if (grepl("Manta", IDLIST, fixed=TRUE) ){
			algorithm = "Manta"
		} else {
			algorithm = "Smoove"
		}

		variant_ID_chr_start_end_type=paste0(chr, "_", AVG_START, "_", AVG_END, "_", type_vcf)

  		# Variant ID, AF_tot, AF_XX, AF_XY, AC_tot, AC_XX, AC_XY, AN_tot, AN_XX, AN_XY, Hom_alt_tot, Hom_alt_XX, Hom_alt_XY
  		# AN_tot : number of 0/0, 0/1 and 1/1 genotypes (avoid counting the ./.)
  		an_total = 2*(sum(GT_table_i == "0/0", na.rm=T) + sum(GT_table_i == "0/1", na.rm=T) + sum(GT_table_i == "1/1", na.rm=T)) 
  		#AC tot
  		ac_total = sum(GT_table_i == "0/1", na.rm=T) + 2*sum(GT_table_i == "1/1", na.rm=T)
  		#AF tot = AC/AN
  		af_total = ac_total/an_total
  		#Number of individus homozygotes for the alternative allele (1/1)
  		hom_alt_total = sum(GT_table_i == "1/1", na.rm=T) 

  		#For XX individuals
  		#For now, make fake false with individuals and sex : 
  		#sex_table =  read.table("sample_sex.tsv", header=TRUE)
  		#Subset the GT_Table for XX individuals
  		XX_Samples = sex_table[sex_table$Sex=="XX",1]
  		XX_GT_table_i = GT_table_i[XX_Samples]
  		# AN_XX
  		an_xx = 2*(sum(XX_GT_table_i == "0/0", na.rm=T) + sum(XX_GT_table_i == "0/1", na.rm=T) + sum(XX_GT_table_i == "1/1", na.rm=T)) 
  		#AC XX
  		ac_xx = sum(XX_GT_table_i == "0/1", na.rm=T) + 2*sum(XX_GT_table_i == "1/1", na.rm=T)
  		#AF X = AC/AN
  		af_xx = ac_xx/an_xx
  		#Number of individus homozygotes for the alternative allele (1/1)
  		hom_alt_xx = sum(XX_GT_table_i == "1/1", na.rm=T) 
  
  		#For XY individuals
  		#Subset the GT_Table for XY individuals
  		XY_Samples = sex_table[sex_table$Sex=="XY",1]
  		XY_GT_table_i = GT_table_i[XY_Samples]
  		# AN_XY
  		an_xy = 2*(sum(XY_GT_table_i == "0/0", na.rm=T) + sum(XY_GT_table_i == "0/1", na.rm=T) + sum(XY_GT_table_i == "1/1", na.rm=T)) 
  		#AC XY
  		ac_xy = sum(XY_GT_table_i == "0/1", na.rm=T) + 2*sum(XY_GT_table_i == "1/1", na.rm=T)
  		#AF X = AC/AN
  		af_xy = ac_xy/an_xy
  		#Number of individus homozygotes for the alternative allele (1/1)
  		hom_alt_xy = sum(XY_GT_table_i == "1/1", na.rm=T) 

		#Some SV are too long and not annotated by the annot workflowm, need to define consequence and gene to avoid error
		if (length(SV_annot_i$Consequence)>1) { 
			#Consequence
			consequence=SV_annot_i$Consequence

                	#gene
			gene=SV_annot_i$SYMBOL
		} else {
			consequence = "NotAnnotated"
			gene = "NotAnnotated"
		}


  		# variant_ID, type, length, chr, pos, ref, alt, cadd_score, cadd_interpr, dbsnp_id, dbsnp_url, UCSC_url, ensembl_url, clinvar_url, gnomad_url
  		# UCSC URL : https://genome.ucsc.edu/cgi-bin/hgTracks?db=<assembly>&highlight=<assembly>.chrM%3A<pos>-<pos>&position=chrM%3A<pos-25>-<pos+25> / hg38 or hg19 db=hg38&highlight=hg38.chrM%3A8602-8602&position=chrM%3A8577-8627
  		# Ensembl_url : https://uswest.ensembl.org/Homo_sapiens/Location/View?r=<chr>%3A<pos-25>-<pos+25> // r=17%3A63992802-64038237
  		# clinvar URL : https://www.ncbi.nlm.nih.gov/clinvar/variation/<VCV>/ // 692920
  		# gnomad_URL : https://gnomad.broadinstitute.org/variant/M-<pos>-<ref>-<alt>?dataset=gnomad_r3 // M-8602-T-C
  

  		# UCSC URL : https://genome.ucsc.edu/cgi-bin/hgTracks?db=<assembly>&highlight=<assembly>.chrM%3A<pos>-<pos>&position=chrM%3A<pos-25>-<pos+25> / hg38 or hg19 db=hg38&highlight=hg38.chrM%3A8602-8602&position=chrM%3A8577-8627
  		ucsc_url=paste0("https://genome.ucsc.edu/cgi-bin/hgTracks?db=",assembly, "&highlight=", assembly, ".chr", chr, "%3A", AVG_START, "-", AVG_END, "&position=chr", chr, "%3A", AVG_START-(0.25*AVG_LEN), "-", AVG_END+(0.25*AVG_LEN))
  
		#Could be added in V2, not part of the current SQL
  		# Ensembl_url : https://uswest.ensembl.org/Homo_sapiens/Location/View?r=<chr>%3A<pos-25>-<pos+25> // r=17%3A63992802-64038237
  		#if (assembly=="GRCh38") {
    		#	ensembl_url=paste0("https://uswest.ensembl.org/Homo_sapiens/Location/View?r=", chr, "%3A", AVG_START-(0.25*AVG_LEN), "-", AVG_END+(0.25*AVG_LEN))
  		#} else if (assembly=="GRCh37") {
    		#	ensembl_url=paste0("https://grch37.ensembl.org/Homo_sapiens/Location/View?r=", chr, "%3A", AVG_START-(0.25*AVG_LEN), "-", AVG_END+(0.25*AVG_LEN))
  		#}  
    
                #Could be added in V2, not part of the current SQL
		# dbsnp_id : From annotation file (SNV_annot_i), "Existing_variation" column
  		#if (grepl("rs", SV_annot_i$Existing_variation)) {
    		#	dbsnp_id = gsub(",.*$", "", SV_annot_i$Existing_variation)
    			# dbsnp_URL : https://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=<rs_number> // rs1556423501
    		#	dbsnp_url=paste0("https://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=", dbsnp_id)
  		#} else {
    		#	dbsnp_id="NA"
    		#	dbsnp_url="NA"
  		#}

		#For SV V1, the gnomAD URL will point to the region view in gnomAD SV
	      	#gnomAD SV is currently only avaialble in GRCh37
		#https://gnomad.broadinstitute.org/region/22-18738600-18747000?dataset=gnomad_sv_r2_1
		gnomad_id = "TBD"
		if (assembly=="GRCh37") {
			gnomad_url_region=paste0("https://gnomad.broadinstitute.org/region/", chr, "-", AVG_START, "-", AVG_END, "?dataset=gnomad_sv_r2_1")
		} else {
                        gnomad_url_region="NA"
                }


		##For complex variations (translocations including a chr2)
		#To do, need to find one in a vcf
		if (type_vcf == "BND") {
			chr2="TBD"
			chr2_pos1="TBD"
			ucsc_url2="TBD"
			gnomad_id2="TBD"
			gnomad_url2="TBD"
			
			#CTX table
                	#variant, chr2, chr2_pos1, ucsc_url2, gnomad_id2, gnomad_url2
                	#TO DO
                	temp_table_ctx_i=cbind(variant, chr2, chr2_pos1, ucsc_url2, gnomad_id2, gnomad_url2)
                	table_ctx=unique(rbind.data.frame(table_ctx, temp_table_ctx_i))
		}

  		### Create tables
  		# SV_IBVL_frequency
  		# Variant ID, AF_tot, AF_XX, AF_XY, AC_tot, AC_XX, AC_XY, AN_tot, AN_XX, AN_XY, Hom_alt_tot, Hom_alt_XX, Hom_alt_XY, qual
  		temp_table_frequ_db_i = cbind(variant, af_total, af_xx, af_xy, ac_total, ac_xx, ac_xy, an_total, an_xx, an_xy, hom_alt_total, hom_alt_xx, hom_alt_xy, quality)
		table_frequ_SV=unique(rbind.data.frame(table_frequ_SV, temp_table_frequ_db_i))

  		# SV_annotation (svs)
  		# variant_ID, chr1, chr1_pos1 (start), chr1_po2 (end), type, length, algorithm, ucsc_url, gnomad_id, gnomad_url
  		temp_table_annot_SV_i = cbind(variant, chr, AVG_START, AVG_END, type_vcf, AVG_LEN, algorithm, ucsc_url, gnomad_id, gnomad_url_region)
		table_annot_SV=unique(rbind.data.frame(table_annot_SV, temp_table_annot_SV_i))


		#sv_consequences table
		#gene, variant, consequence (intronic, intergenic, etc)
		#No transcript associated to SV, so no ensembl or Refseq
		# If there is several consequences on the same line (separrated by a coma), create one line per consequence
		temp_table_sv_consequence_i = cbind(gene, variant, consequence)
		show(i)
		show(SV_annot_i)
		show(length(SV_annot_i$Consequence))
		show(temp_table_sv_consequence_i)
		table_sv_consequence=unique(rbind.data.frame(table_sv_consequence, temp_table_sv_consequence_i))
		table_sv_consequence_split=separate_rows(table_sv_consequence, consequence, sep = ",")

	}       
       
	### Write table slots
    	write.table(table_frequ_SV, file=paste0("table_frequ_SV_slot", j), quote=FALSE, row.names = FALSE)
    	write.table(table_annot_SV, file=paste0("table_annot_SV_slot", j), quote=FALSE, row.names = FALSE)
	write.table(table_sv_consequence_split, file=paste0("table_SV_consequence_slot", j), quote=FALSE, row.names = FALSE)
	if (nrow(table_ctx)>0){
		write.table(table_ctx, file=paste0("table_ctx_slot", j), quote=FALSE, row.names = FALSE)

	}
}

### Write tables 
# SV_IBVL_frequency
list_frequ_tables_slots <- list.files(pattern = paste0("table_frequ_SV_slot"))
tables_frequ_slots=lapply(list_frequ_tables_slots, read.table, header=TRUE)
combined_tables_frequ_slots=do.call(rbind, tables_frequ_slots)
write.table(combined_tables_frequ_slots, file=paste0("genomic_ibvl_frequencies_", var_type, "_", chromosome, ".tsv"), quote=FALSE, row.names = FALSE, sep="\t")
file.remove(list_frequ_tables_slots)

# SV_annotation
list_annot_tables_slots <- list.files(pattern = paste0("table_annot_SV_slot"))
tables_annot_slots=lapply(list_annot_tables_slots, read.table, header=TRUE)
combined_tables_annot_slots=do.call(rbind, tables_annot_slots)
colnames(combined_tables_annot_slots)=c("variant", "chr1", "chr1_pos1", "chr1_pos2", "type", "length", "algorithm", "ucsc_url", "gnomad_id", "gnomad_url")
write.table(combined_tables_annot_slots, file=paste0("svs_", var_type, "_", chromosome,".tsv"), quote=FALSE, row.names = FALSE, sep="\t")
file.remove(list_annot_tables_slots)


# Variants_consequences table
list_variant_consequence_tables_slots <- list.files(pattern = paste0("table_SV_consequence_slot"))
tables_variant_consequence_slots=lapply(list_variant_consequence_tables_slots, read.table, header=TRUE)
combined_tables_variant_consequence_slots=do.call(rbind, tables_variant_consequence_slots)
combined_tables_variant_consequence_slots_unique=unique(combined_tables_variant_consequence_slots)
write.table(combined_tables_variant_consequence_slots_unique, file=paste0("sv_consequences_", var_type, "_", chromosome,".tsv"), quote=FALSE, row.names = FALSE, sep="\t")
file.remove(list_variant_consequence_tables_slots)

# ctx table
list_ctx_tables_slots <- list.files(pattern = paste0("table_ctx_slot"))
tables_ctx_slots=lapply(list_ctx_tables_slots, read.table, header=TRUE)
combined_tables_ctx_slots=do.call(rbind, tables_ctx_slots)
if (length(nrow(combined_tables_ctx_slots))>0) {
	write.table(combined_tables_ctx_slots, file=paste0("svs_ctx_", chromosome,".tsv"), quote=FALSE, row.names = FALSE, sep="\t")
}
file.remove(list_ctx_tables_slots)

#Variants
variants_table=cbind(as.data.frame(unique(SV_raw_annotaton_file[,c("Uploaded_variation")])), "SV")
colnames(variants_table)=c("variant_id", "var_type")
write.table(variants_table, file=paste0("variants_sv_", var_type, "_", chromosome,".tsv"), quote=FALSE, row.names = FALSE, sep="\t")

# Gene table
# Short name, (NCBI gene Id, may not be necessary anymore)
gene_table=as.data.frame(unique(SV_raw_annotaton_file[,c("SYMBOL")]))
colnames(gene_table)=c("short_name")
gene_table=as.data.frame(gene_table[!grepl("^-$", gene_table$short_name),])
colnames(gene_table)=c("short_name")
write.table(gene_table, file=paste0("genes_sv_", var_type, "_", chromosome, ".tsv"), quote=FALSE, row.names = FALSE, sep="\t")


