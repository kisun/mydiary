#!/bin/bash -l
#SBATCH -J feelnc-filter
#SBATCH -o out_feelnc-filter_%j.txt
#SBATCH -e err_feelnc-filter_%j.txt
#SBATCH -t 4:00:00
#SBATCH -n 4
#SBATCH --mem-per-cpu=8000

gtf="/wrk/kipokh/DONOTREMOVE/Data/Genomes/oar3.0/oar31_87.gtf"
dir="/wrk/kipokh/DONOTREMOVE/Data/projects/mRNA/cl-end/results/hisat2/"
genome="/wrk/kipokh/DONOTREMOVE/Data/Genomes/oar3.0/oar31.fa"
outdir="/wrk/kipokh/DONOTREMOVE/Data/projects/mRNA/cl-end/results/feelnc/"

module load bioconda/3
source activate bioconda_tools



FEELnc_filter.pl \
-i ${dir}/cl-end_merged.gtf \
-a $gtf \
-b transcript_biotype=protein_coding \
> ${outdir}/candidate_lncRNA.gtf

