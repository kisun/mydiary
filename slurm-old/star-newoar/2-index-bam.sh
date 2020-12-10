#!/bin/bash -l
#SBATCH -J saidx
#SBATCH -o out_saidx_%A_%a.txt
#SBATCH -e err_saidx_%A_%a.txt
#SBATCH -t 1:00:00
#SBATCH -n 1
#SBATCH --array=1-42
#SBATCH --mem-per-cpu=4000
#SBATCH --cpus-per-task=4

cd /wrk/kipokh/DONOTREMOVE/Data/projects/mRNA/cl-end/results/star_rambo

bam=$(ls *.bam | sed -n "$SLURM_ARRAY_TASK_ID"p)
echo "Analysing sample file: $bam"
ls -l $bam

samtools index $bam 

