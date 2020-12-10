#!/bin/bash -l
#SBATCH -J gffcompare
#SBATCH -o out_gffcompare_%j.txt
#SBATCH -e err_gffcompare_%j.txt
#SBATCH -t 01:00:00
#SBATCH -n 1
#SBATCH --mem-per-cpu=24000

gtf="/wrk/kipokh/DONOTREMOVE/Data/Genomes/oar3.0/oar31_87.gtf"
dir="/wrk/kipokh/DONOTREMOVE/Data/projects/mRNA/cl-end/results/hisat2/"

module load bioconda/3
source activate bioconda_tools

cd $dir
gffcompare -r $gtf -o cl-end_gffcompare cl-end_merged.gtf
