#!/bin/bash

#-------------------------------------------------------
#SBATCH -J setup_ref_genome_multi
#SBATCH -N 1 
#SBATCH -n 1
#SBATCH -p development
#SBATCH -o logs/setup_ref_genome.o
#SBATCH -e logs/setup_ref_genome.e
#SBATCH -t 01:00:00
#-------------------------------------------------------
ml tacc-singularity
# Script for setting up genome files for bismark
# This version combines the Amil and Asel genomes
module unload xalt

mkdir genomes_for_bismark
cd genomes_for_bismark


REF_DIR=$STOCKYARD/tagmap-share/genomes
REMOTE_AMIL=$REF_DIR/Amil_v2.1/GCF_013753865.1_Amil_v2.1_genomic.fna
REMOTE_AMIL_GFF=$REF_DIR/Amil_v2.1/GCF_013753865.1_Amil_v2.1_genomic.gff
REMOTE_ASEL=$REF_DIR/Asel_v1.0/asel.fasta
REMOTE_ASEL_GFF=$REF_DIR/Asel_v1.0/asel.gff

# Where to save this

# Directory with prep scripts
SCRIPTS=$PWD/scripts
HERE=$PWD
#### WOrkhorse function

function make_genome {
  local NAME=$1
  local GENOME="${NAME}_lambda.fasta"
  local GFF="${NAME}.gff"
  local CORES=${2:-64}
  mkdir $NAME
  cd $NAME
  mv ../$GENOME .
  # Append the lambda genome
  # Modify some lambda tags and attach to Amil
  sed 's/>NC_001416.1/>lambda_phage/' ../GCF_000840245.1_ViralProj14204_genomic.fna >> $GENOME
  
  # run the bismark prep script
  bismark_genome_preparation --parallel $CORES .
  
  # Add the GFF
  mv ../$GFF .
  
    #generate gene windows
  $SCRIPTS/gff_to_bed4.py -gff $GFF -feature gene -IDstring ID -o geneBoundaries.bed &
  
  #generate exon windows
  $SCRIPTS/gff_to_bed4.py -gff $GFF -feature CDS -IDstring ID -o cdsBoundaries.bed &
  
  #generate promoter windows
  $SCRIPTS/gff_to_promoter_bed.py -gff $GFF -bp 1000 -IDstring ID -o promoterBoundaries.bed &
  
  #generate for tss windows
  $SCRIPTS/gff_to_tssWindow_bed.py -gff $GFF -bp 250 -IDstring ID -o tssBoundaries.bed & 

  biopython $SCRIPTS/fasta_sequence_characters.py -fa $GENOME > chrLengths.txt 
  
  #generate 1Kb windows
  bedtools makewindows -g chrLengths.txt -w 1000 | awk 'BEGIN{OFS="\t"}{print $0,$1"-"$2"-"$3}' > windowBoundaries_1kb.bed &
  
  #generate 500 bp windows
  bedtools makewindows -g chrLengths.txt -w 500 | awk 'BEGIN{OFS="\t"}{print $0,$1"-"$2"-"$3}' > windowBoundaries_500bp.bed 
  cd $HERE
}

# Generate the windows; this requires more cores so should be run serially
function make_window {
  local NAME=$1
  local GENOME="${NAME}_lambda.fasta"
  local GFF="${NAME}.gff"
  local CORES=${2:-128}
  cd $NAME
  r-tidy-opt $SCRIPTS/windows_around_genes_from_gff.R --cores $CORES --i $GFF --o ${GFF}_around_genes_1kb.bed
  cd $HERE
}

function make_reference {
  local NAME=$1
  local GENOME="${NAME}_lambda.fasta"
  local GFF="${NAME}.gff"
  local GENOME_FOLDER="${PWD}/${NAME}"
  local OUT_FILE="${SCRATCH}/hybrid_methylation/genomes/${NAME}.sh"
#------- MAKE CALLABLE REFERENCE FILE -------#
#we end up needing lots of files for the reference genome,
#so it's easy to assmeble all the paths in a single .sh file,
#so you can switch to that reference easily from $SCRTACH/
#(this is helpful when working with multiple species at once)


echo "\
GENOME_FOLDER=${GENOME_FOLDER}
GENOME_PATH=${GENOME_FOLDER}/${GENOME}
GFF_PATH=${GENOME_FOLDER}/${GFF}
GENE_ID=ID
exonWindowFile=${GENOME_FOLDER}/cdsBoundaries.bed
geneWindowFile=${GENOME_FOLDER}/geneBoundaries.bed
promoterWindowFile=${GENOME_FOLDER}/promoterBoundaries.bed
tssWindowFile=${GENOME_FOLDER}/tssBoundaries.bed
window1KbFile=${GENOME_FOLDER}/windowBoundaries_1kb.bed
window500bpFile=${GENOME_FOLDER}/windowBoundaries_500bp.bed
aroundGeneWindowFile=${GENOME_FOLDER}/${GFF}_around_genes_1kb.bed" > $OUT_FILE

}


# We want the lambda genome tacked on to test conversion efficiency
# I'm not really sure why this is
# Ask someone about it?

echo "Copying the reference genomes from remote source"
cp $REMOTE_AMIL Amil_lambda.fasta
cp $REMOTE_ASEL Asel_lambda.fasta
cp $REMOTE_AMIL_GFF Amil.gff
cp $REMOTE_ASEL_GFF Asel.gff

# create Amilsel
cat Amil_lambda.fasta Asel_lambda.fasta > Amilsel_lambda.fasta
cat Amil.gff Asel.gff > Amilsel.gff

echo "Downloading the Lambda Genome"
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/840/245/GCF_000840245.1_ViralProj14204/GCF_000840245.1_ViralProj14204_genomic.fna.gz
gunzip GCF_*_ViralProj14204_genomic.fna.gz

echo "Creating Bismark versions of genomes"

make_genome Amil 32 &
make_genome Asel 32 &
make_genome Amilsel 32 &

wait

echo "Creating windows"

# Make the windows
make_window Amil
make_window Asel
make_window Amilsel

mkdir $SCRATCH/hybrid_methylation/genomes

make_reference Amil 
make_reference Asel &
make_reference Amilsel &

wait
