#!/bin/bash -l
#SBATCH -J s2sb
#SBATCH -o out_s2sb_%A_%a.txt
#SBATCH -e err_s2sb_%A_%a.txt
#SBATCH -t 2:00:00
#SBATCH --array=1-42
#SBATCH --mem-per-cpu=24000
#SBATCH --cpus-per-task=2

cd /wrk/kipokh/DONOTREMOVE/Data/projects/mRNA/cl-end/results/hisat2/bam/

file=$(ls *.sam | sed -n "$SLURM_ARRAY_TASK_ID"p )
base_name="${file%%.*}"
#module load biokit

#RUN samtools

samtools view -Su $file | samtools sort - -o ${base_name}_sorted.bam
