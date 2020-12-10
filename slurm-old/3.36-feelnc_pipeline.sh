#!/bin/bash -l
#SBATCH -J feelnc
#SBATCH -o out_feelnc_%j.txt
#SBATCH -e err_feelnc_%j.txt
#SBATCH -t 12:00:00
#SBATCH -n 4
#SBATCH --mem-per-cpu=8000

gtf="/wrk/kipokh/DONOTREMOVE/Data/Genomes/oar3.0/oar31_87.gtf"
dir="/wrk/kipokh/DONOTREMOVE/Data/projects/mRNA/cl-end/results/hisat2/"
genome="/wrk/kipokh/DONOTREMOVE/Data/Genomes/oar3.0/oar31.fa"
outdir="/wrk/kipokh/DONOTREMOVE/Data/projects/mRNA/cl-end/results/feelnc/"

module load bioconda/3
source activate bioconda_tools



FEELnc_pipeline.sh \
--candidate=${dir}/cl-end_merged.gtf \
--reference=$gtf \
--genome=$genome \
--outname=cl-end \
--outdir=$outdir

