#!/bin/bash -l
#SBATCH -J kallisto
#SBATCH -o out_kallisto_%A_%a.txt
#SBATCH -e err_kallisto_%A_%a.txt
#SBATCH -t 1:00:00
#SBATCH -n 4
#SBATCH --array=1-42
#SBATCH --mem-per-cpu=4000
#SBATCH --cpus-per-task=3

cd /wrk/kipokh/DONOTREMOVE/Data/projects/mRNA/cl-end/results/trimgalore

R1=$(ls *_R1_val_1.fq.gz | sed -n "$SLURM_ARRAY_TASK_ID"p )
R2=$(echo ${R1} | sed 's/_R1_val_1.fq.gz/_R2_val_2.fq.gz/')
base_name="${R1%%_R1*}"


#RUN kallisto
kallisto \
quant \
--index=/wrk/kipokh/DONOTREMOVE/Data/Genomes/oar3.0/kallisto-index \
--output-dir=/wrk/kipokh/DONOTREMOVE/Data/projects/mRNA/cl-end/results/kallisto/${base_name} \
--threads=4 \
--plaintext \
$R1 $R2
