#!/bin/bash -l
#SBATCH -J stargenome
#SBATCH -o out_stargenome_%A_%a.txt
#SBATCH -e err_stargenome_%A_%a.txt
#SBATCH -t 02:00:00
#SBATCH -n 4
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=12000
#
#Load biopython-env

cd /wrk/kipokh/DONOTREMOVE/Data/Genomes/oar_rambo_v1.0/


#######################################################
## run star genome ##
#######################################################

STAR \
--runMode genomeGenerate \
--genomeDir /wrk/kipokh/DONOTREMOVE/Data/Genomes/oar_rambo_v1.0/star-genome \
--genomeFastaFiles /wrk/kipokh/DONOTREMOVE/Data/Genomes/oar_rambo_v1.0/oar_rambo.fa \
--runThreadN 4 \
--sjdbGTFfile /wrk/kipokh/DONOTREMOVE/Data/Genomes/oar_rambo_v1.0/oar_rambo.gtf \
--sjdbOverhang 99

