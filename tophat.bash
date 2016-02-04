#!/bin/bash -l
#SBATCH -J TopHat
#SBATCH -o output_%j.txt
#SBATCH -e errors_%j.txt
#SBATCH -t 12:00:00
#SBATCH -n 1
#SBATCH --array=1-2
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2000

module load biokit
cd /wrk/kipokh/Data/mRNA/mRaw/mRaw2/two/

#set up file names for a subjob
i=$(ls *_R1_001.fastq.gz | sed -n "$SLURM_ARRAY_TASK_ID"p )
b=$(echo ${i} | sed 's/_R1_001.fastq.gz/_R2_001.fastq.gz/')
dir="${i%%_*}"

#Run tophat
tophat -p 8 -o /wrk/kipokh/Data/mRNA/mAlign/Sheep_Sheep/mAlign2S/ENS75/$dir \
-G /wrk/kipokh/Data/Genomes/Sheep/ENS75/Annotation/Oarv31Ens75.gtf \
/wrk/kipokh/Data/Genomes/Sheep/ENS75/Sequence/Oarv31Ens75 $i $b
