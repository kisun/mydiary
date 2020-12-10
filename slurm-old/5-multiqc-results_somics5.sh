#!/bin/bash -l
#SBATCH -J multiqc
#SBATCH -o out_multiqc_%j.txt
#SBATCH -e err_multiqc_%j.txt
#SBATCH -t 07:00:00
#SBATCH -n 1
#SBATCH --mem-per-cpu=24000

cd /wrk/kipokh/DONOTREMOVE/Data/projects/mRNA/cl-end/results/

module load bioconda/3
source activate bioconda_tools

multiqc -d fastqc trimgalore hisat2 salmon rseqc-hisat2 -n multiqc_somics5
