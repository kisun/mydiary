#!/bin/bash -l
#SBATCH --job-name=star-genome
#SBATCH --output=inembp_out.txt
#SBATCH --error=inembp_err.txt
#SBATCH --account=project_2003269
#SBATCH --partition=small
#SBATCH --time=05:00:00
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=6000

module load biokit

STAR \
--runThreadN 4 \
--runMode genomeGenerate \
--genomeDir /scratch/project_2003269/genome/cow/fasta/starGenome_btaARSUCD1.2 \
--genomeFastaFiles /scratch/project_2003269/genome/cow/fasta/bta-ARS_UCD1.2/index/bta-ARS_UCD1.2_toplevel.fa \
--sjdbGTFfile /scratch/project_2003269/genome/cow/gtf/bta_ARS-UCD1.2.gtf \
--sjdbOverhang 149





