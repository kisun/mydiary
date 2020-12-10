#!/bin/bash -l
#SBATCH -J salmon
#SBATCH -o out_salmon_%A_%a.txt
#SBATCH -e err_salmon_%A_%a.txt
#SBATCH -t 4:00:00
#SBATCH -n 4
#SBATCH --array=1-42
#SBATCH --mem-per-cpu=4000
#SBATCH --cpus-per-task=4

cd /wrk/kipokh/DONOTREMOVE/Data/RawData/mRNA/III/cl-end

R1=$(ls *_R1_001.fastq.gz | sed -n "$SLURM_ARRAY_TASK_ID"p )
R2=$(echo ${R1} | sed 's/_R1_001.fastq.gz/_R2_001.fastq.gz/')
base_name="${R1%%_*}"

#module load biokit

#RUN salmon
salmon \
quant \
-i /wrk/kipokh/DONOTREMOVE/Data/Genomes/oar3.0/salmon-index \
-l A \
-1 $R1 \
-2 $R2 \
-o /wrk/kipokh/DONOTREMOVE/Data/projects/mRNA/cl-end/results/salmon/${base_name}_quant


