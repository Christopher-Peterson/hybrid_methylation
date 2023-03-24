#!/bin/bash
# Helper script for SNPsplit_genome_preparation

VCF=$1 # input VCF file, filtered
GENOME_NAME=${2:-Amil}
source genome.sh # GENOME_FOLDER
PARENTS=`basename "${VCF/_filtered.vcf/}"`

# Set up a subdirectory to work in
ROOTDIR=$PWD
WD=masked_genomes/$PARENTS
mkdir $WD # If this direc exists, the run will fail.
ln -s $ROOTDIR/$VCF $WD/parents.vcf
ln -s $ROOTDIR/scripts $WD/scripts

cd $WD

SCRIPT=scripts/SNPsplit_genome_preparation

# List the strains in the file; this will be a reverse-order array
STRAINS=( `$SCRIPT --vcf_file parents.vcf --list_strains 2> /dev/null` )

# Run the scripts
$SCRIPT --vcf_file parents.vcf \
	  --strain "${STRAINS[1]}" --strain2 "${STRAINS[0]}" \
	    --reference_genome $GENOME_FOLDER --dual_hybrid \
	      --nmasking --genome_build $GENOME_NAME

# Should I aslo add Bismark Prep to this?

# Enter the masked genome file
cd A?_A?_dual*masked

# Remove 'chr' leader from all resulting filenames
for i in chr*.fa; do
	   newName="${i/chr/}"
	      mv "$i" "$newName"
      done



      # RUn genome prep script
      bismark_genome_preparation --parallel 20 .


      # make a shell file

       cd ..
       mkdir -p ../registry # In the masked_genomes
       GENOME_FOLDER=`ls $PWD/A?_A?_dual*masked`

      # Copy genome to $WORK
#      mkdir -p $WORK/hybrid_methylation/masked_genomes/$PARENTS
#      cp -r Bisulfite_Genome/ $WORK/hybrid_methylation/masked_genomes/$PARENTS & 
      # FIgure out offspring from parents?

      mkdir -p $ROOTDIR/masked_genomes/registry

      # Write out the
      echo "GENOME_FOLDER=$WORK/hybrid_methylation/masked_genomes/$PARENTS" >  ${ROOTDIR}/masked_genomes/registry/${PARENTS}.sh

