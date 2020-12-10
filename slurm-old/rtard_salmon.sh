#!/bin/bash -l
#SBATCH -J salmon
#SBATCH -o out_salmon_%A_%a.txt
#SBATCH -e err_salmon_%A_%a.txt
#SBATCH -t 4:00:00
#SBATCH -n 4
#SBATCH --array=1-24
#SBATCH --mem-per-cpu=4000
#SBATCH --cpus-per-task=4

cd /wrk/kipokh/DONOTREMOVE/Data/RawData/mRNA/ArcArk_RNAseq/rta/merged/

R1=$(ls *_R1.fastq.gz | sed -n "$SLURM_ARRAY_TASK_ID"p )
R2=$(echo ${R1} | sed 's/_R1.fastq.gz/_R2.fastq.gz/')
base_name="${R1%%_R1*}"

#module load biokit

#RUN salmon

/homeappl/home/kipokh/appl_taito/mytools/Salmon-0.8.2_linux_x86_64/bin/salmon \
quant \
-i /wrk/kipokh/DONOTREMOVE/Data/Genomes/rta1.0/rtard/regulation/rta_salmon_index \
-l A \
-1 $R1 \
-2 $R2 \
-o /wrk/kipokh/DONOTREMOVE/Data/mRNA/rtard/salmon/salmon-counts/${base_name}_quant


