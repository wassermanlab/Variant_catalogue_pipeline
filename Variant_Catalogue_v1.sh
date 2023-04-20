#!/bin/sh

## CPU Usage
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4G
#SBATCH --time=200:30:00
#SBATCH -p silent_q
#SBATCH --mail-user=melsiddieg@cmmt.ubc.ca
# #SBATCH --mail-type=FAIL

## Output and Stderr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.error

source /mnt/common/SILENT/Act3/conda/miniconda3/etc/profile.d/conda.sh
Nextflow=/mnt/common/Precision/NextFlow/nextflow
module load singularity
prof=$1
$Nextflow run main.nf -profile $prof -resume -with-trace -with-report -with-timeline  -with-dag flowchart.png
