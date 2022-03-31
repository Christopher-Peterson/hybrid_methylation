#!/bin/bash

#-------------------------------------------------------
#SBATCH -J setup_ref_genome
#SBATCH -N 1 
#SBATCH -n 1
#SBATCH -p development
#SBATCH -o logs/setup_ref_genome.o
#SBATCH -e logs/setup_ref_genome.e
#SBATCH -t 01:00:00
#-------------------------------------------------------

# Script for setting up genome files for bismark

REF_DIR=$STOCKYARD/tagmap-share/genomes
REMOTE_GENOME=$REF_DIR/Amil_v2.1/GCF_013753865.1_Amil_v2.1_genomic.fna
REMOTE_GFF=$REF_DIR/Amil_v2.1/GCF_013753865.1_Amil_v2.1_genomic.gff
# Local files; to be made
GFF=Amil.coding.gff
GENOME=Amil_lambda.fasta

# Where to save this
GENOME_FOLDER=$REF_DIR/Amil_v2.1_ForBismark
# Shell script file to save references to
OUT_FILE=$REF_DIR/amilleporaForBismark.sh

# Directory with prep scripts
SCRIPTS=`pwd`/scripts
# Number of cores for the windowing R script
WINDOW_CORES=128 

ml tacc-singularity
module unload xalt
##### Setup the genome file
mkdir genomes_for_bismark

# We want the lambda genome tacked on to test conversion efficiency
# I'm not really sure why this is
# Ask someone about it?

cd genomes_for_bismark
echo "Copying the reference genome from remote source"
cp $REMOTE_GENOME $GENOME

echo "Downloading the Lambda Genome"
mkdir tmp
cd tmp
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/840/245/GCF_000840245.1_ViralProj14204/GCF_000840245.1_ViralProj14204_genomic.fna.gz
gunzip GCF_*_ViralProj14204_genomic.fna.gz
cd .. # genomes for bismark

# Modify some lambda tags and attach to Amil
sed 's/>NC_001416.1/>lambda_phage/' tmp/GCF_000840245.1_ViralProj14204_genomic.fna >> $GENOME

# Clean up
rm tmp/*
rmdir tmp

# make sure only the Amil_lambda is in there
FILES_IN_FOLDER=$(ls | wc -l)
echo "Files in the bismark genome folder: $FILES_IN_FOLDER (there should only be one:)"

# run the bismark prep script
bismark_genome_preparation --parallel 64 . #genomes_for_bismark 

FILES_IN_FOLDER=$(ls Bisulfite_Genome/CT_conversion| wc -l)
echo "Files in the CT conversion directory: $FILES_IN_FOLDER (there should be some):"

# Create windows in the genome ###########################
echo "Creating specific genomic windows"
cp $REMOTE_GFF $GFF # download GFF file
#generate gene windows
$SCRIPTS/gff_to_bed4.py -gff $GFF -feature gene -IDstring ID -o geneBoundaries.bed &

#generate exon windows
$SCRIPTS/gff_to_bed4.py -gff $GFF -feature CDS -IDstring ID -o cdsBoundaries.bed &

#generate promoter windows
$SCRIPTS/gff_to_promoter_bed.py -gff $GFF -bp 1000 -IDstring ID -o promoterBoundaries.bed &

#generate for tss windows
$SCRIPTS/gff_to_tssWindow_bed.py -gff $GFF -bp 250 -IDstring ID -o tssBoundaries.bed &


# BIOPY_SIF=/work2/04386/crpeters/ls6/singularity/biopython_latest.sif
# TIDYV_SIF=/work2/04386/crpeters/ls6/singularity/r-tidyverse-optparse_4.1.2.sif

# chmod +x $SCRIPTS/fasta_sequence_characters.py
# Jump into the container to run this
biopython $SCRIPTS/fasta_sequence_characters.py -fa $GENOME > chrLengths.txt 

#generate 1Kb windows
bedtools makewindows -g chrLengths.txt -w 1000 | awk 'BEGIN{OFS="\t"}{print $0,$1"-"$2"-"$3}' > windowBoundaries_1kb.bed &

#generate 500 bp windows
bedtools makewindows -g chrLengths.txt -w 500 | awk 'BEGIN{OFS="\t"}{print $0,$1"-"$2"-"$3}' > windowBoundaries_500bp.bed &

wait
#1Kb windows around genes
# singularity shell -H /home $TIDYV_SIF
# GFF=Amil.coding.gff3
# GENOME=Amil_lambda.fasta

# Add N_CORES option
r-tidy-opt $SCRIPTS/windows_around_genes_from_gff.R --cores $WINDOW_CORES --i $GFF --o ${GFF}_around_genes_1kb.bed
# The initial version of this was painfully slow; I've re-written it if I need to repeat it


#------- MAKE CALLABLE REFERENCE FILE -------#
#we end up needing lots of files for the reference genome,
#so it's easy to assmeble all the paths in a single .sh file,
#so you can switch to that reference easily from $SCRTACH/
#(this is helpful when working with multiple species at once)

echo "Copying genome to $GENOME_FOLDER"
mkdir $GENOME_FOLDER
cp -r *  $GENOME_FOLDER

echo "\
GENOME_FOLDER=${GENOME_FOLDER}
GENOME_PATH=${GENOME_FOLDER}/Amil_lambda.fasta
GFF_PATH=${GENOME_FOLDER}/Amil.coding.gff3
GENE_ID=ID
exonWindowFile=${GENOME_FOLDER}/cdsBoundaries.bed
geneWindowFile=${GENOME_FOLDER}/geneBoundaries.bed
promoterWindowFile=${GENOME_FOLDER}/promoterBoundaries.bed
tssWindowFile=${GENOME_FOLDER}/tssBoundaries.bed
window1KbFile=${GENOME_FOLDER}/windowBoundaries_1kb.bed
window500bpFile=${GENOME_FOLDER}/windowBoundaries_500bp.bed
aroundGeneWindowFile=${GENOME_FOLDER}/Amil.coding.gff3_1kb_around_genes.bed" > $OUT_FILE

# Create a symlink to hybrid_methylation/genome.sh
cds hybrid_methylation
ln -s $OUT_FILE genome.sh
