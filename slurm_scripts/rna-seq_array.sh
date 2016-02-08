#!/bin/bash -l
#SBATCH -J array_job
#SBATCH -o array_job_out_%A_%a.txt
#SBATCH -e array_job_err_%A_%a.txt
#SBATCH -t 00:20:00
#SBATCH --mem-per-cpu=4000
#SBATCH --array=1-2
#SBATCH -n 1
#SBATCH -p test

# Load modules
# module load biokit
# module load R
# module load biopython-env

# Make directories for result files

mkdir -p results-fastqc
mkdir -p results-trimmomatic
mkdir -p results-htseq
mkdir -p results-rseqc
mkdir -p results-bam
mkdir -p results-tophat

# Choose the filename from list
sample_file=$(sed -n "$SLURM_ARRAY_TASK_ID"p fastq.list)

# Remove the extension to use the base name as a name for folders etc. Note that the base name (sample name) should not contain dots.
base_name="${sample_file%%.*}"

# FastQC
fastqc -o results-fastqc ${sample_file}

# Trimmomatic
trimmomatic SE -threads 1 -phred33 ${sample_file} results-trimmomatic/${base_name}_trimmed.fq.gz TRAILING:5 MINLEN:50

# TopHat2
tophat2 -o results-tophat/${base_name} --transcriptome-index=tophat-indexes/transcriptome/hg38chr19 tophat-indexes/hg38chr19 results-trimmomatic/${base_name}_trimmed.fq.gz
mv results-tophat/${base_name}/accepted_hits.bam results-bam/${base_name}.bam

# samtools
samtools index results-bam/${base_name}.bam

# RseQC
bam_stat.py -i results-bam/${base_name}.bam 2> results-rseqc/${base_name}-stats.txt
read_distribution.py -r refseq_19.bed -i results-bam/${base_name}.bam > results-rseqc/${base_name}-distribution.txt
geneBody_coverage.py -r refseq_19.bed -i results-bam/${base_name}.bam -o results-rseqc/${base_name}
junction_annotation.py -r refseq_19.bed -i results-bam/${base_name}.bam -o results-rseqc/${base_name} 2> results-rseqc/${base_name}-junction-counts.txt
junction_saturation.py -r refseq_19.bed -i results-bam/${base_name}.bam -o results-rseqc/${base_name}

# HTSeq
htseq-count -f bam --stranded=no results-bam/${base_name}.bam tophat-indexes/hg38chr19.gtf | head -n -5 > results-htseq/${base_name}-counts.tsv
