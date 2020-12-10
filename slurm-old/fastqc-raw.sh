#!/bin/bash -l
#SBATCH -J fastqc_raw
#SBATCH -o out_qc-raw_%A_%a.txt
#SBATCH -e err_qc-raw_%A_%a.txt
#SBATCH -t 02:00:00
#SBATCH -n 1
#SBATCH --array=1-84
#SBATCH --mem-per-cpu=4000

#module load biokit: run this command before executing this script. 

#####go to the files directory

cd /wrk/kipokh/DONOTREMOVE/Data/RawData/mRNA/III/cl-end/
sample_file=$(ls *.fastq.gz | sed -n "$SLURM_ARRAY_TASK_ID"p)

echo "Analysing sample file: $sample_file"
ls -l $sample_file

#######################################################
## Fastqc ##
#######################################################

fastqc -o /wrk/kipokh/DONOTREMOVE/Data/projects/mRNA/cl-end/results/fastqc ${sample_file} 

#mv *.txt err-out/fastqc/
