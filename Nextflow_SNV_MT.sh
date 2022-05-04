#!/bin/sh

## CPU Usage
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4G
#SBATCH --time=200:30:00
#SBATCH -p silent_q

## Output and Stderr
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.error

module load singularity
nextflow run Nextflow_SNV_MT_211220.nf -profile GRCh37 -resume -with-trace -with-report -with-timeline  -with-dag flowchart.png
