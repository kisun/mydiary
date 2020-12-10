#!/bin/bash -l
#SBATCH -J fc-star
#SBATCH -o out_fc-star_%A_%a.txt
#SBATCH -e err_fc-star_%A_%a.txt
#SBATCH -t 02:00:00
#SBATCH -n 1
#SBATCH --array=1-42
#SBATCH --mem-per-cpu=16000



gtf="/wrk/kipokh/DONOTREMOVE/Data/Genomes/oar3.0/oar31_87.gtf"
dir="/wrk/kipokh/DONOTREMOVE/Data/projects/mRNA/cl-end/results/star/"

cd $dir
bam=$(ls *Aligned.sorted*.bam | sed -n "$SLURM_ARRAY_TASK_ID"p)
bamid="${bam%%Aligned.sorted*}"
echo "Analysing sample file: $bam"
ls -l $bam

featureCounts -a $gtf -g gene_id -o ${bamid}_gene.featureCounts.txt -p $bam
featureCounts -a $gtf -g gene_biotype -o ${bamid}_biotype.featureCounts.txt -p $bam
cut -f 1,7 ${bamid}_biotype.featureCounts.txt > ${bamid}_biotype_counts.txt
