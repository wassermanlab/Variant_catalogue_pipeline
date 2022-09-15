#!/bin/bash

########################################################################
### FILE        install.sh
###
### DESCRIPTION This is the install.sh script for the CAFE pipeline project.
###		It will use singularity images, and conda environments, for tooling
###
### Phillip Richmond
### 2022-09-14
########################################################################

# set -e

########
# Usage:
# bash install.sh -d </path/to/database/install/directory/>
# DESC This code block gets the database directory from the input options and stores it as DB_DIR
#########
usage()
{
    echo "usage: bash install.sh -d </path/to/database/install/directory/>"
    echo "example: $ bash install.sh -d /mnt/common/DATABASES/"
}


# Get directory for install
while getopts d: flag
do
    case "${flag}" in
        d) DB_DIR=${OPTARG};;
    esac
done

if [ -z "$DB_DIR" ]
then
    usage
    exit
fi

echo "Databases will be installed in: $DB_DIR";

# Get directory that this pipeline lives in
CAFE_Pipeline_Dir=$PWD




############
# nextflow #
############

##########################################
# DESC This function deploys nextflow executable
# ARGS This function doesnt require any arguments
# RSLT Downloads nextflow executable from internet and changes to executable permission
##########################################

function buildNextflow()
{
	# make directory where we'll deploy nextflow
	mkdir -p $CAFE_Pipeline_Dir/Nextflow
	cd $CAFE_Pipeline_Dir/Nextflow
	
	# standard fetch from interweb and run
	curl -s https://get.nextflow.io | bash

	# change permissions
	chmod ugo=rwx ./nextflow

	# check
	./nextflow --help
	cd $CAFE_Pipeline_Dir
}
buildNextflow



###################
# buildMiniConda3 #
###################

##########################################
# DESC This function deploys miniconda3
# ARGS This function doesnt require any arguments
# RSLT Downloads miniconda3 installer and deploys a miniconda3 directory inside this repo
##########################################

# This function fetches and builds miniconda 
function buildMiniConda3() 
{
	Miniconda_Install_Dir=$CAFE_Pipeline_Dir/miniconda3
	
	# get installer script
        wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
        bash Miniconda3-latest-Linux-x86_64.sh -b -p $Miniconda_Install_Dir

	# add conda to your path
        source ${Miniconda_Install_Dir}/etc/profile.d/conda.sh
        cd $CAFE_Pipeline_Dir

	# Check the install
	conda --help

	# Install libmamba solver
	conda install -y conda-libmamba-solver=22.8

	# cleanup
	rm Miniconda3-latest-Linux-x86_64.sh
}
buildMiniConda3


##################
# buildCondaCAFE #
##################


##########################################
# DESC This function uses miniconda3 to deploy a single conda environment 
# ARGS This function doesnt require any arguments
# RSLT Will create an env which can be activated with conda activate CAFE_pipeline
##########################################
# This will build an environment for each of the necessary tools being run from conda, including
# bcftools=1.9
# htslib=1.9
# bwa=0.7.17
# vcftools
# haplocheck
# bwa
# see the yml for the rest
function buildCAFE() 
{
	Miniconda_Install_Dir=$CAFE_Pipeline_Dir/miniconda3
	
	# add conda to path
	source ${Miniconda_Install_Dir}/etc/profile.d/conda.sh

	# install CAFE pipeline tools into an env
	conda env create --experimental-solver=libmamba -f $CAFE_Pipeline_Dir/CAFE_pipeline_conda.yml 
	
	# test it
	cd $CAFE_Pipeline_Dir
	conda activate CAFE_pipeline

	# install hail
	conda deactivate
	conda create --experimental-solver=libmamba -n hail -c bioconda -y hail=0.2.61
	conda activate hail

}
buildCAFE


###############
# Pull Images #
###############

##########################################
# DESC This function pull images from online repositories using wget or singularity pull, 
# 	and puts them into a dir called singularity/
# 	This function relies on singularity being installed, and the ability to run: `singularity pull`
# 	This function pulls images needed for CAFE including:
# 		- DeepVariant
# 		- FastQC
# 		- GATK
# 		- GLNexus
# 		- SV Tools from BrentP: Jasmine, Paragraph, 
# 		- Manta
# 		- MosDepth
# 		- multiQC
# 		- Smoove
# 		- VEP
# ARGS This function doesnt require any arguments
# RSLT Will create an env which can be activated with conda activate CAFE_pipeline
##########################################


function PullImages()
{
	Dir=$PWD
	
	# make a directory to deploy the images
	Singularity_Dir=$Dir/singularity
	mkdir -p $Singularity_Dir
	cd $Singularity_Dir

	# DeepVariant
	singularity pull deepvariant-1.2.0.sif docker://google/deepvariant:1.2.0
	
	# FastQC
	wget -cq -O fastqc-0.11.9.sif https://depot.galaxyproject.org/singularity/fastqc%3A0.11.9--hdfd78af_1

	# GATK
	wget -cq -O gatk-4.2.0.sif https://depot.galaxyproject.org/singularity/gatk4%3A4.2.0.0--0

	# GLnexus
	wget -cq -O glnexus-1.4.1.sif https://depot.galaxyproject.org/singularity/glnexus%3A1.4.1--h40d77a6_0

	# SV tools from brentp
	singularity pull rare_disease_sv-0.1.2.sif docker://brentp/rare-disease-sv:v0.1.2

	# Manta
	wget -cq -O manta-1.6.sif https://depot.galaxyproject.org/singularity/manta%3A1.6.0--py27_0

	# Mosdepth
	wget -cq -O mosdepth-0.3.2.sif https://depot.galaxyproject.org/singularity/mosdepth%3A0.3.2--h01d7912_0

	# multiQC
	wget -cq -O multiqc-1.9.sif https://depot.galaxyproject.org/singularity/multiqc%3A1.9--py_1

	# R
	singularity pull bioconductor_docker-3.15.sif docker://bioconductor/bioconductor_docker:RELEASE_3_15

	# Smoove
	wget -cq -O smoove-0.2.8.sif https://depot.galaxyproject.org/singularity/smoove%3A0.2.8--h9ee0642_0

	# VEP
	wget -cq -O vep-92.4.sif https://depot.galaxyproject.org/singularity/ensembl-vep%3A92.4--htslib1.7_0
	
}
PullImages


exit

## Below is in development for fetching databases

function createDatabaseDirectory()
{
	DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
	OPT_DIR=${DIR}/opt
	cd $DIR
	
	# We need bwa and samtools, so load the environment we created above
	MINI_CONDA_INSTALL_DIR=$OPT_DIR/miniconda3
	source ${MINI_CONDA_INSTALL_DIR}/etc/profile.d/conda.sh
	conda activate opt/AnnotateVariantsEnvironment

	mkdir -p $DB_DIR

	# Get Reference Genome Files 
	## GRCh37
	mkdir -p $DB_DIR/GRCh37/
	cd $DB_DIR/GRCh37/
	wget http://www.bcgsc.ca/downloads/genomes/9606/hg19/1000genomes/bwa_ind/genome/GRCh37-lite.fa
	wget http://www.bcgsc.ca/downloads/genomes/9606/hg19/1000genomes/bwa_ind/genome/GRCh37-lite.fa.fai

	bwa index GRCh37-lite.fa

	## GRCh38
	mkdir -p $DB_DIR/GRCh38/
	wget ftp://ftp.ensembl.org/pub/release-96/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
	gunzip -c Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz > GRCh38-lite.fa

	samtools faidx GRCh38-lite.fa
	bwa index GRCh38-lite.fa

}


function getDatabases()
{
	echo "not yet"
}

