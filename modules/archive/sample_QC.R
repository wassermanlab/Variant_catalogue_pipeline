#R script
#Created by Solenne Correard in December 2021
#Owned by the Silent Genomes Project Activity 3 team
#Developped to build the IBVL, a background variant library

#Script to define the sex of each individual and create graph to define the thresholds for sample QC
#Sex is defined on three things in gnomAD, same here
#The current pipeline uses a rough F-stat cutoff of 0.5 to split samples into the XX and XY categories. The final X and Y ploidy cutoffs are then determined from the means and standard deviations of those XX and XY distributions. Sex was assigned based on the following cutoffs:
#In gnomAD, sample QC is done based on some hard filters. The threasholds may difffer from gnomAD to the IBVL, so graph representing the distribution of the parameters (Numbers of SNPs per sample, sumber of singletons per sample, het/hom ratio) are created to see if gnomAD thresholds are relevant or if other threasholds should be implemented in the IBVL.


library(ggplot2)
library('ggplot2')

#Sex definition
#XY:
#normalized X coverage < 1.29 &
#normalized Y coverage > 0.1 &
#normalized Y coverage < 1.16
#XX:
#normalized X coverage > 1.45 &
#normalized X coverage < 2.4 &
#normalized Y coverage < 0.1

args <- commandArgs(trailingOnly = TRUE)

assembly = (args[1])
batch = (args[2])
run = (args[3])

#Load plink file with F values (used to define the sex of each sample)
plink_F_file =read.table(args[4], header=TRUE)
#plink_F_file =read.table("work/35/96a2ce6cf4ebbed3596f786a1928a5/DeepVariant_GLnexus_Run_20211220_plink_sex.sexcheck", header=TRUE)

#Load the file with number of singletons per sample created using vcftools
#singleton_file = read.table(args[5], header=TRUE)
#singleton_file = read.table("work/b4/b9555c434b5aac68b983d69646c42c/singleton_per_ind.tsv", header=TRUE)

#Load the file with number of SNPs per sample created using GATK
#SNPs_file = read.table(args[6], sep=":")
#SNPs_file = read.table("/mnt/scratch/SILENT/Act3/Processed/Individual/GRCh37/Batch_DryRun/Run_20211220/QC/Count_variants/SNP_count.tsv", sep=":")

#Load the bcftools stat file. It includes also the number of singletons per sample (Number differ from number given by vcftools, need to compare both) and the het/hom ratio
#bcftools_stats= read.table(args[7])
#bcftools_stats= read.table("work/bd/df922cd588aafa56ef445e82931848/bcftools_stats.tsv")
#colnames(bcftools_stats)=c("PSC", "[2]id", "sample", "[4]nRefHom", "nNonRefHom", "nHets", "[7]nTransitions", "[8]nTransversions", "[9]nIndels", "[10]average_depth", "nSingletons", "[12]nHapRef", "[13]nHapAlt", "[14]nMissing")


table_QC=data.frame()
#For each sample
for (i in 1: nrow(plink_F_file)) {
	sample = plink_F_file$IID[i]
	plink_F = plink_F_file$F[i]	
	mosdepth_file_i = read.table(paste0(sample, "_sorted.mosdepth.summary.txt"), header= TRUE)
#        mosdepth_file_i = read.table(paste0("/mnt/scratch/SILENT/Act3/Processed/Individual/GRCh37/Batch_DryRun/Run_20211220/QC/",sample,"_sorted/Mosdepth/",sample, "_sorted.mosdepth.summary.txt"), header= TRUE)
	X_coverage = mosdepth_file_i$mean[mosdepth_file_i$chrom=="X"]
	Y_coverage = mosdepth_file_i$mean[mosdepth_file_i$chrom=="Y"]
	chr20_coverage = mosdepth_file_i$mean[mosdepth_file_i$chrom=="20"]
	mean_coverage = mosdepth_file_i$mean[mosdepth_file_i$chrom=="total"]

	X_normalized_coverage = 2*(X_coverage / chr20_coverage)
	Y_normalized_coverage = Y_coverage / chr20_coverage

	# gnomad : The current pipeline uses a rough F-stat cutoff of 0.5 to split samples into the XX and XY categories.
	# plink : By default, F estimates smaller than 0.2 yield female calls, and values larger than 0.8 yield male calls. 
	if (plink_F<0.2 & X_normalized_coverage > 1.45 & X_normalized_coverage < 2.4 & Y_normalized_coverage < 0.1) {
		Sex = "XX"
	} else if (plink_F>0.8 & X_normalized_coverage < 1.29 & Y_normalized_coverage > 0.1 & Y_normalized_coverage < 1.16) {
		Sex = "XY"
	} else {
		Sex = "ambiguous"
	}

#	singleton = singleton_file$X1[singleton_file$INDV==sample]
#	SNPs = SNPs_file$V2[SNPs_file$V1==paste0(sample, "_sorted_gatk_count_variants")]
#	singleton_bcftools = bcftools_stats$nSingletons[bcftools_stats$sample==sample]
#	het_hom_ratio = bcftools_stats$nHets[bcftools_stats$sample==sample] / bcftools_stats$nNonRefHom[bcftools_stats$sample==sample]

#	temp_table_QC = cbind(sample, Sex, mean_coverage, singleton, SNPs, singleton_bcftools, het_hom_ratio, plink_F, X_normalized_coverage, Y_normalized_coverage)
        temp_table_QC = cbind(sample, Sex, mean_coverage, plink_F, X_normalized_coverage, Y_normalized_coverage)
	table_QC=rbind.data.frame(table_QC, temp_table_QC)
}

write.table(table_QC, file="QC_sample.tsv", quote=FALSE, row.names = FALSE, sep="\t")

##Create graph
#Sex inference graph : x : chr X relative ploidy (0.5 to 3), Y : Chr Y relative ploidy (0 to 1.6)
# Add a rectangle in the bakground to identify variant that fall into the XX or XY area. Variant with ambiguous sex will be outside of the rectangles
# XX rectangle : 0 < Y < 0.1, 1.45 < X < 2.4
# XY rectangle : 0.1 < Y < 1.16 , 0 < X < 1.29

sex_graph = ggplot(table_QC) + 
  geom_point(aes(as.numeric(X_normalized_coverage), as.numeric(Y_normalized_coverage))) + 
  ggtitle ("Sex determination") +
  xlab ("chr X normalized coverage") +
  ylab ("chr Y normalized coverage") +
  geom_rect(data=data.frame(xmin = 1.45, xmax = 2.4, ymin = -Inf, ymax = 0.1),
            aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="pink", alpha=0.5) +
  geom_rect(data=data.frame(xmin = -Inf, xmax = 1.29, ymin = 0.1, ymax = 1.16),
            aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="blue", alpha=0.5)
ggsave("sex_graph.pdf")


#Singleton graph
# Distribution of the number of singleton by individual
# gnomAD :outliers were filtered using the following cutoffs: Number of singletons: > 100k
# gnomAD threashold: 100,000 

#Number of singletons counted using vcftools
#pdf(file="singleton_graph.pdf", width = 4, height = 4)
#hist(as.numeric(table_QC$singleton), breaks=12, xlab="Number of singleton SNV per sample",
#     xlim = c(50000,3000000), main="Singleton distribution by vcftools")
#  abline(v = 100000, col="blue")
#dev.off()

#Number of singletons counted using bcftools
#pdf(file="singleton_bcftools_graph.pdf", width = 4, height = 4)
#hist(as.numeric(table_QC$singleton_bcftools), breaks=12, xlab="Number of singleton SNV per sample",
#     main="Singleton distribution by bcftools")
#  abline(v = 100000, col="blue")
#dev.off()

#Coverage graph
# Distribution of the coverage by individual (on the whole genome, while gnomAD does it only based on chr20)
# gnomAD : samples were filtered if they had a mean coverage on chromosome 20 < 15X
# gnomAD threashold : 15

#pdf(file="coverage_graph.pdf", width = 4, height = 4)
#hist(as.numeric(table_QC$mean_coverage), breaks=12, xlab="Mean coverage per sample",
#     xlim = c(0,50), main = "Mean coverage distribution")
#abline(v = 15, col="blue")  	
#dev.off()

# Total number of SNP graph
# Distribution of the total number of SNPs by individuals
# gnomAD : Outliers were filtered using the following cutoffs: Number of SNVs: < 2.4M or > 3.75M
# gnomAD threasholds at < 2.4M and > 3.75M

#pdf(file="SNPs_count_graph.pdf", width = 4, height = 4)
#hist(as.numeric(table_QC$SNPs), breaks=12, xlab="Total number of SNPs per sample",
#       xlim = c(200000,10000000), main = "SNPs count distribution")
#  abline(v = 2400000, col="blue")+
#    abline(v = 3750000, col="blue")
#dev.off()


# Ratio number of heterozygous variants / number of homozygous variants
# A good ratio would be around 2, a ratio way higher could indicate contamination by a sample from the same specie
# gnomAD : Outliers were filtered using the following cutoffs: Ratio of heterozygous to homozygous variants: > 3.3
# gnomAD threashold at 3.3

#pdf(file="het_hom_ratio_graph.pdf", width = 4, height = 4)
#hist(as.numeric(table_QC$het_hom_ratio), breaks=12, xlab="Het/Hom ratio per sample",
#       xlim = c(0,5), main = "Het/Hom distribution")
#  abline(v = 3.3, col="blue")
#dev.off()


