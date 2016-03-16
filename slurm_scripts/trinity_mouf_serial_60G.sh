#!/bin/bash -l
#SBATCH -J Mouf_Trin
#SBATCH -o output_%j.txt
#SBATCH -e errors_%j.txt
#SBATCH -t 3-00:00:00
#SBATCH -n 1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16000
#

module load biokit

#######################################################
##  Run Trinity to Generate Transcriptome Assemblies ##
#######################################################
## use jellyfish
Trinity --seqType fq --max_memory 60G --trimmomatic --normalize_reads --left \
/wrk/kipokh/Data/Mouflon/all_left.fastq.gz --right /wrk/kipokh/Data/Mouflon/all_right.fastq.gz \
--SS_lib_type RF --CPU $SLURM_CPUS_PER_TASK --output /wrk/kipokh/Data/Mouflon/trinity_mouf_results_serial_60G \
--grid_conf $TRINITY_HOME/hpc_conf/taito.slurm


