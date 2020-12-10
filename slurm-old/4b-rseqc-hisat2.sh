#!/bin/bash -l
#SBATCH -J rseqc-hs2
#SBATCH -o out_rseqc-hs2_%A_%a.txt
#SBATCH -e err_rseqc-hs2_%A_%a.txt
#SBATCH -t 08:00:00
#SBATCH -n 1
#SBATCH --array=1-42
#SBATCH --mem-per-cpu=24000

bed12="/wrk/kipokh/DONOTREMOVE/Data/Genomes/oar3.0/oar31_87.bed"
gtf="/wrk/kipokh/DONOTREMOVE/Data/Genomes/oar3.0/oar31_87.gtf"
outdir="/wrk/kipokh/DONOTREMOVE/Data/projects/mRNA/cl-end/results/rseqc-hisat2"
dir="/wrk/kipokh/DONOTREMOVE/Data/projects/mRNA/cl-end/results/hisat2/bam/"

cd $dir
bam=$(ls *_sorted.bam | sed -n "$SLURM_ARRAY_TASK_ID"p)
bamid="${bam%%_sorted*}"
echo "Analysing sample file: $bam"
ls -l $bam

samtools index $bam 
infer_experiment.py -i $bam -r $bed12 > $outdir/$bamid.infer_experiment.txt
RPKM_saturation.py -i $bam -r $bed12 -d $strandRule -o $outdir/$bamid.RPKM_saturation
junction_annotation.py -i $bam -o $outdir/$bamid.rseqc -r $bed12
bam_stat.py -i $bam > $outdir/$bamid.bam_stat.txt
junction_saturation.py -i $bam -o $outdir/$bamid.rseqc -r $bed12 
inner_distance.py -i $bam -o $outdir/$bamid.rseqc -r $bed12
geneBody_coverage.py -i $bam -o $outdir/$bamid.rseqc -r $bed12
read_distribution.py -i $bam -r $bed12 > $outdir/$bamid.read_distribution.txt
read_duplication.py -i $bam -o $outdir/$bamid.read_duplications
