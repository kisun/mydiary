#!/bin/bash -l
#SBATCH -J merge-fc-hs2
#SBATCH -o out_merge-fc-hs2_%A_%a.txt
#SBATCH -e err_merge-fc-hs2_%A_%a.txt
#SBATCH -t 01:00:00
#SBATCH -n 1
#SBATCH --mem-per-cpu=16000

module load python-env

dir="/wrk/kipokh/DONOTREMOVE/Data/projects/mRNA/cl-end/results/hisat2/"

./merge_featurecounts.py -o $dir/merged_gene_counts.txt -i $dir/*_gene.featureCounts.txt
