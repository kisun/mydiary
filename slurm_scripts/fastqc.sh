#!/bin/bash -l
#SBATCH -J fastqc
#SBATCH -o output_%j.txt
#SBATCH -e errors_%j.txt
#SBATCH -t 03:00:00
#SBATCH -n 1
#SBATCH --nodes=1
#SBATCH --array=1-43
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2000
#

module load biokit

#######################################################
## Fastqc ##
#######################################################

fastqc *.fastq.gz -o /wrk/kipokh/Data/mRNA/mQuali/MQC_3/CL
