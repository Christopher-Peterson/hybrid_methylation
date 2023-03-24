#!/bin/bash

#SBATCH -J index_ref
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p normal
#SBATCH -o logs/index_ref.o
#SBATCH -e logs/index_ref.e
#SBATCH -t 0:45:00

  
module unload xalt
cd ref_genome
ln -Lfs Amil_lambda.fasta Amil_lambda.fa
BatMeth2 index -g Amil_lambda.fa