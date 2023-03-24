#!/bin/bash

#-------------------------------------------------------
#SBATCH -J index_masked_genomes  # sed_* will be replaced w/ sed on creation of experiment
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p normal
#SBATCH -o logs/index_masked_genomes.o
#SBATCH -e logs/index_masked_genomes.e
#SBATCH -t 1:20:00  

#------------------------------------------------------

# Filter the bcf files to only include variants w/ two samples

# ml launcher
module unload xalt
source genome.sh
cd masked_genomes
mkdir -p indexed



source genome.sh
CHR_ORDER=( `cat $GENOME_FOLDER/chrLengths.txt | grep -Po "^.+\t" ` )

function concat_masked_genome {
  local parent=${1}
  local GENOME_FILES=""
  for chr in ${CHR_ORDER[@]}; do 
    local GENOME_FILES="$GENOME_FILES ${parent}/*hybrid*masked/${chr}.N-masked.fa"
  done
  #echo "Genome files: $GENOME_FILES"
  cat $GENOME_FILES > indexed/${parent}.fa
}


for dir in A??; do
  concat_masked_genome $dir
done

cd indexed
for genome in A??.fa; do
  bowtie2-build -f --threads 20 $genome ${genome/.fa/} &
done




# ln -lfs Amil_lambda.fasta Amil_lambda.fa
# BatMeth2 index -g Amil_lambda.fa

# BatMeth2 index -g masked_genomes/A13/full_N-masked.fa

# Get the order of the chromosomes
# CHR_ORDER=( `cat $GENOME_FOLDER/chrLengths.txt | grep -Po "^.+\t" ` )
# 
# function concat_masked {
#   local dir=${1}
#   local GENOME_FILES=""
#   for chr in ${CHR_ORDER[@]}; do 
#     local GENOME_FILES="$GENOME_FILES ${dir}/*hybrid*masked/${chr}.N-masked.fa "
#   done
#   echo "Genome files: $GENOME_FILES"
#   cat $GENOME_FILES > indexed/${dir}.fa
# }
# for dir in A??; do
#   concat_masked $dir &
# done
# wait

# Now they're concated, index them
cd indexed
for genome in A.fa; do
  BatMeth2 index -g $genome &
done
wait
# 
# # Setup the job file
# # export LAUNCHER_JOB_FILE=jobs/mask_genomes
# 
# # Assemble the job file
# > $LAUNCHER_JOB_FILE
# for vcf in vcf/concat/A??_filtered.vcf; do
#   echo "scripts/run_SNPsplit_genome_prep.sh $vcf $GENOME_NAME" >> $LAUNCHER_JOB_FILE
# done
# 
# export LAUNCHER_WORKDIR=`pwd`
# $LAUNCHER_DIR/paramrun
# 
