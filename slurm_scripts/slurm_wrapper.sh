#!/bin/bash
############################################################################
## Description:
## The wrapper script to run the RNA-Seq workflow step by step using SLURM
##
## Author: Kisun Pokharel (kisun.pokharel@helsinki.fi
##
############################################################################
#!/bin/bash
##note first few steps may be run only once for one project, comment and 
##change their dependencies accordingly for related jobs
Genome=$(sbatch getgenome.sh| cut -f 4 -d' ')
echo $Genome
Genes=$(sbatch getgenes.sh| cut -f 4 -d' ')
echo $Genes
Bowtie2Build=$(sbatch bowtie2build.sh -d afterok:$Genome:$Genes| cut -f 4 -d' ') 
echo $Bowtie2Build
Fastqc1=$(sbatch fastqc_pretrim.sh --array=1-30| cut -f 4 -d' ') #change the --array value based on number of samples
echo $Fastqc1
Trimmomatic=$(sbatch -d afterok:$Fastqc1 trimmomatic.sh --array=1-30| cut -f 4 -d' ')
echo $Trimmomatic
Fastqc2=$(sbatch -d afterok:$Trimmomatic fastqc_posttrim.sh --array=1-30| cut -f 4 -d' ')
echo $Fastqc2
Tophat2=$(sbatch -d afterok:$Fastqc2:$Bowtie2Build tophat2.sh --array=1-30| cut -f 4 -d' ')
echo $Tophat2
Cufflinks=$(sbatch -d afterok:$Tophat2 cufflinks.sh --array=1-30| cut -f 4 -d' ')
echo $Cufflinks
Htseqcount=$(sbatch -d afterok:$Tophat2 htseqcount.sh --array=1-30| cut -f 4 -d' ')
echo $Htseqcount
exit 0

#and so on
