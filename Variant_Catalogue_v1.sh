#!/bin/sh

## CPU Usage
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4G
#SBATCH -p defq
#SBATCH --mail-user=prichmond@bcchr.ca
#SBATCH --mail-type=ALL

## Output and Stderr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.error

source /mnt/common/SILENT/Act3/conda/miniconda3/etc/profile.d/conda.sh
Nextflow=/mnt/common/Precision/NextFlow/nextflow
module load singularity

$Nextflow run Variant_Catalogue_v1.nf -profile GRCh38 -resume -with-trace -with-report -with-timeline  -with-dag flowchart.png
