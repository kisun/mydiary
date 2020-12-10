#!/bin/bash -l
#SBATCH -J stringtie
#SBATCH -o out_stie_%A_%a.txt
#SBATCH -e err_stie_%A_%a.txt
#SBATCH -t 2:00:00
#SBATCH --array=1-42
#SBATCH --mem-per-cpu=4000
#SBATCH --cpus-per-task=4

cd /wrk/kipokh/DONOTREMOVE/Data/projects/mRNA/cl-end/results/hisat2/bam/

file=$(ls *.bam | sed -n "$SLURM_ARRAY_TASK_ID"p )
base_name="${file%%.*}"
#module load biokit

#RUN stringtie

stringtie ${file} \
-p $SLURM_CPUS_PER_TASK \
-G /wrk/kipokh/DONOTREMOVE/Data/Genomes/oar3.0/oar31_87.gtf \
-l ${base_name} \
-e \
-A \
-B \
--merge 
