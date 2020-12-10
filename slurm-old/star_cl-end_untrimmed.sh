#!/bin/bash -l
#SBATCH -J star
#SBATCH -o out_star_%A_%a.txt
#SBATCH -e err_star_%A_%a.txt
#SBATCH -t 6:00:00
#SBATCH -n 4
#SBATCH --array=1-42
#SBATCH --mem-per-cpu=4000
#SBATCH --cpus-per-task=4

cd /wrk/kipokh/DONOTREMOVE/Data/RawData/mRNA/III/cl-end

R1=$(ls *_R1_001.fastq.gz | sed -n "$SLURM_ARRAY_TASK_ID"p )
R2=$(echo ${R1} | sed 's/_R1_001.fastq.gz/_R2_001.fastq.gz/')
base_name="${R1%%_.*}"

STAR \
--genomeDir /wrk/kipokh/DONOTREMOVE/Data/Genomes/oar3.0/star-genome1 \
--readFilesIn $R1 $R2 \
--readFilesCommand zcat \
--runThreadN 4 \
--sjdbGTFfile /wrk/kipokh/DONOTREMOVE/Data/Genomes/oar3.0/oar31_87.gtf \
--outFileNamePrefix /wrk/kipokh/DONOTREMOVE/Data/projects/mRNA/cl-end/results/star/${base_name} \
--quantMode GeneCounts \
--twopassMode Basic

