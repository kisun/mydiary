#!/bin/bash -l
#SBATCH -J fc-hs2
#SBATCH -o out_fc-hs2_%A_%a.txt
#SBATCH -e err_fc-hs2_%A_%a.txt
#SBATCH -t 02:00:00
#SBATCH -n 1
#SBATCH --array=1-42
#SBATCH --mem-per-cpu=16000



gtf="/wrk/kipokh/DONOTREMOVE/Data/Genomes/oar3.0/oar31_87.gtf"
dir="/wrk/kipokh/DONOTREMOVE/Data/projects/mRNA/cl-end/results/hisat2/bam/"

cd $dir
bam=$(ls *_sorted.bam | sed -n "$SLURM_ARRAY_TASK_ID"p)
bamid="${bam%%_sorted*}"
echo "Analysing sample file: $bam"
ls -l $bam

featureCounts -a $gtf -g gene_id -o ${bamid}_gene.featureCounts.txt -p $bam
featureCounts -a $gtf -g gene_biotype -o ${bamid}_biotype.featureCounts.txt -p $bam
cut -f 1,7 ${bamid}_biotype.featureCounts.txt > ${bamid}_biotype_counts.txt
