#!/bin/sh

## CPU Usage
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4G
#SBATCH --time=200:30:00
#SBATCH -p defq

## Output and Stderr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.error

source /conda/miniconda3/etc/profile.d/conda.sh
Nextflow=n/NextFlow/nextflow
module load singularity

$Nextflow run main.nf -profile GRCh37 -resume -with-trace -with-report -with-timeline  -with-dag flowchart.png
