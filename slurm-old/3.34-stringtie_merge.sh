#!/bin/bash -l
#SBATCH -J stringtie-merge
#SBATCH -o out_stringtie-merge_%j.txt
#SBATCH -e err_stringtie-merge_%j.txt
#SBATCH -t 02:00:00
#SBATCH -n 1
#SBATCH --mem-per-cpu=24000

gtf="/wrk/kipokh/DONOTREMOVE/Data/Genomes/oar3.0/oar31_87.gtf"
dir="/wrk/kipokh/DONOTREMOVE/Data/projects/mRNA/cl-end/results/hisat2/"

module load bioconda/3
source activate bioconda_tools

cd $dir
stringtie --merge -G $gtf -o cl-end_merged.gtf cl-end_gtf_list.txt
