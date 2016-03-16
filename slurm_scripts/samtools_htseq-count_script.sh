#!/bin/bash -l
# created: Feb 13, 2014 11:06 AM
# author: kipokh
#SBATCH -J samtools
#SBATCH -o out
#SBATCH -e error
#SBATCH -n 1
#SBATCH -t 04:00:00
#SBATCH --mem-per-cpu=8000
#SBATCH --mail-type=END
#SBATCH --mail-user=kisun.pokharel@helsinki.fi

# commands to manage the batch script
#   submission command
#     sbatch [script-file]
#   status command
#     squeue -u kipokh
#   termination command
#     scancel [jobid]

# For more information
#   man sbatch
#   more examples in Taito guide in http://datakeskus.csc.fi/web/guest/taito-user-guide

# copy this script to your terminal and then add your commands here

#example run commands

#srun ./my_serial_program

