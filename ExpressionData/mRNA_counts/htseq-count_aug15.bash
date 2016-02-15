#!/bin/bash -l
# created: Feb 13, 2014 11:06 AM
# author: kipokh
#SBATCH -J H-C
#SBATCH -o out_%j.txt
#SBATCH -e error_%j.txt
#SBATCH -n 1
#SBATCH --array=1-39
#SBATCH --cpus-per-task=1
#SBATCH -t 08:00:00
#SBATCH --mem-per-cpu=24000

##To use HTSeq-Count, the tool needs to be manually installed in application directory. 
##But before, installing, the python environment should be loaded with the command below:

module load biopython-env

cd  /wrk/kipokh/Data/mRNA/mAlign/Sheep_Sheep/mAlign2S/ENS75/Aug15/allaccepted
i=$(ls *_nsort.sam | sed -n "$SLURM_ARRAY_TASK_ID"p)
b=$(echo ${i} | sed 's/_nsort.sam/_counts/')
#samtools sort -n $i $b
#samtools view -h $b > $c
htseq-count -s no $i /wrk/kipokh/Data/Genomes/Sheep/ENS75/Annotation/Oarv31Ens81.gtf > /wrk/kipokh/Data/mRNA/mExpre/Sheep_Sheep/mExpre2S/Aug15Ens81/$b
