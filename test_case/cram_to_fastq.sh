#Script in order to generate fastq from the cram for the 1000 genomes dataset
#Details from the cram format and compression process by the IGSR are available here : http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/README_using_1000genomes_cram.md

Sample=$(basename $filename | cut -f1 -d ".")
INCRAM=${Sample}.final.cram
INBAM=${Sample}.final.bam
INBAM_SORTED=${Sample}.namesorted.bam
OUTFQ1=${Sample}.R1.fastq.gz
OUTFQ2=${Sample}.R2.fastq.gz

# Samtools view
samtools view -hbu -T GRCh38_full_analysis_set_plus_decoy_hla.fa $INCRAM | \
samtools sort -n -o $INBAM_SORTED 

# Convert to fastq
samtools fastq \
		-1 $OUTFQ1 -2 $OUTFQ2 \
		-0 /dev/null \
		-s /dev/null \
		-n $INBAM_SORTED

# Clean
rm $INBAM_SORTED
