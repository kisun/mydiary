#!/bin/bash -l
#SBATCH -J stringtie
#SBATCH -o out_stringtie_%A_%a.txt
#SBATCH -e err_stringtie_%A_%a.txt
#SBATCH -t 03:00:00
#SBATCH --array=1-42
#SBATCH --mem-per-cpu=4000
#SBATCH --cpus-per-task=4



gtf="/wrk/kipokh/DONOTREMOVE/Data/Genomes/oar3.0/oar31_87.gtf"
dir="/wrk/kipokh/DONOTREMOVE/Data/projects/mRNA/cl-end/results/hisat2/"

cd $dir
bam=$(ls *_sorted.bam | sed -n "$SLURM_ARRAY_TASK_ID"p)
bamid="${bam%%_sorted*}"
echo "Analysing sample file: $bam"
ls -l $bam


stringtie $bam \
-p $SLURM_CPUS_PER_TASK \
-o ${bamid}_transcripts.gtf \
-l ${bamid} \
-v \
-G $gtf \
-A \
-e \
-b \
-C
