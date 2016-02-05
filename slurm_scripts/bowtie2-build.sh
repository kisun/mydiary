#!/bin/bash -l
#SBATCH -J Bowtie-build
#SBATCH -o output_%j.txt
#SBATCH -e errors_%j.txt
#SBATCH -t 03:00:00
#SBATCH -n 1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2000
#

module load biokit

#######################################################
## Bowtie build ##
#######################################################

bowtie-build RefGenome.fa RefGenome

