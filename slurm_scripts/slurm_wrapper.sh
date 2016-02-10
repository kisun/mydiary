#!/bin/bash
############################################################################
## Description:
## The wrapper script to run the RNA-Seq workflow
##
## Author: Kisun Pokharel (kisun.pokharel@helsinki.fi
##
############################################################################
#!/bin/bash 

Fastqc1=$(sbatch fastqc_pretrim.sh | cut -f 4 -d' ')
echo $Fastqc1
Trimmomatic=$(sbatch -d afterok:$Fastqc1 trimmomatic.sh | cut -f 4 -d' ')
echo $Trimmomatic
Fastqc2=$(sbatch -d afterok:$Trimmomatic fastqc_posttrim.sh | cut -f 4 -d' ')
echo $Fastqc2
Tophat2=$(sbatch -d afterok:$Fastqc2 tophat2.sh | cut -f 4 -d' ')
echo $Tophat2
Cufflinks=$(sbatch -d afterok:$Tophat2 cufflinks.sh | cut -f 4 -d' ')
echo $Cufflinks
Htseqcount=$(sbatch -d afterok:$Tophat2 htseqcount.sh | cut -f 4 -d' ')
echo $Htseqcount
exit 0

#and so on
