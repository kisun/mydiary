#!/bin/bash -l
#SBATCH -J trimgalore
#SBATCH -o out_trimgalore_%A_%a.txt
#SBATCH -e err_trimgalore_%A_%a.txt
#SBATCH -t 4:00:00
#SBATCH -n 2
#SBATCH --array=1-42
#SBATCH --mem-per-cpu=4000
#SBATCH --cpus-per-task=4



cd /wrk/kipokh/DONOTREMOVE/Data/RawData/mRNA/III/cl-end-ren

R1=$(ls *_R1.fastq.gz | sed -n "$SLURM_ARRAY_TASK_ID"p )
R2=$(echo ${R1} | sed 's/_R1.fastq.gz/_R2.fastq.gz/')



## Trimgalore ##
#######################################################

trim_galore --paired --fastqc --illumina --gzip -o /wrk/kipokh/DONOTREMOVE/Data/projects/mRNA/cl-end/results/trimgalore $R1 $R2
