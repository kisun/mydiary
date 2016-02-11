#!/bin/bash -l
#SBATCH -J TopHat
#SBATCH -o output_%j.txt
#SBATCH -e errors_%j.txt
#SBATCH -t 12:00:00
#SBATCH -n 1
#SBATCH --array=1-2
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2000

wget -P ReferenceGenome
wget -P ReferenceAnnotation
