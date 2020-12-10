#!/bin/bash -l
#SBATCH -J star
#SBATCH -o out_star_%A_%a.txt
#SBATCH -e err_star_%A_%a.txt
#SBATCH -t 5:00:00
#SBATCH -n 4
#SBATCH --array=1-42
#SBATCH --mem-per-cpu=6000
#SBATCH --cpus-per-task=4

cd /wrk/kipokh/DONOTREMOVE/Data/projects/mRNA/cl-end/results/trimgalore

R1=$(ls *_R1_val_1.fq.gz | sed -n "$SLURM_ARRAY_TASK_ID"p )
R2=$(echo ${R1} | sed 's/_R1_val_1.fq.gz/_R2_val_2.fq.gz/')
base_name="${R1%%_R1*}"

STAR \
--genomeDir /wrk/kipokh/DONOTREMOVE/Data/Genomes/oar_rambo_v1.0/star-genome \
--readFilesIn $R1 $R2 \
--readFilesCommand zcat \
--runThreadN 4 \
--sjdbGTFfile /wrk/kipokh/DONOTREMOVE/Data/Genomes/oar_rambo_v1.0/oar_rambo.gtf \
--outFileNamePrefix /wrk/kipokh/DONOTREMOVE/Data/projects/mRNA/cl-end/results/star_rambo/${base_name} \
--outSAMtype BAM SortedByCoordinate \
--quantMode GeneCounts \
--twopassMode Basic


