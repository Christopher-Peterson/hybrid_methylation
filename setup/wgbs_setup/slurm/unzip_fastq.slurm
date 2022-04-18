#!/bin/bash

#-------------------------------------------------------
#SBATCH -J fastq_unzip
#SBATCH -N 1 
#SBATCH -n 12
#SBATCH -p development
#SBATCH -o logs/unzip.o
#SBATCH -e logs/unzip.e
#SBATCH -t 00:30:00

#------------------------------------------------------
cds hybrid_methylation/wgbs_setup/raw_data/fastqs

# Unzip with pigz (faster than regular gunzip)
for file in lane*.gz
do pigz -d $file &
done

wait
