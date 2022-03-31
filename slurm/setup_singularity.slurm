#!/bin/bash

#-------------------------------------------------------
#SBATCH -J setup_singularity
#SBATCH -N 1 
#SBATCH -n 1
#SBATCH -p development
#SBATCH -o logs/setup_ref_genome.o
#SBATCH -e logs/setup_ref_genome.e
#SBATCH -t 00:20:00
#-------------------------------------------------------

# This installs singularity images on TACC
# You need singularity_mask.sh (from the docker folder) installed in $WORK/singularity and for $WORK/singularity/bin to be added to $PATH in .bashrc

cdw singularity
# Install singularity images
ml tacc-singularity
singularity pull -F docker://crpeters/trim_galore:0.6.7 &
singularity pull -F docker://crpeters/bedtools &
singularity pull -F docker://biopython/biopython:latest &
singularity pull -F docker://crpeters/r-tidyverse-optparse:4.1.2 &
# singularity pull -F docker://crpeters/seqmonk:latest &

# add in methylseq
# Install nextcore
curl -fsSL get.nextflow.io | bash
mv nextflow ~/bin/
# Install methylseq 
nextflow run nf-core/methylseq -profile test,singularity 

# Create the binary wrappers for each of them
source singularity_mask.sh

singularity_mask bedtools $WORK/singularity/bedtools_latest.sif
singularity_mask bismark $WORK/nextflow/work/singularity/nfcore-methylseq-1.6.1.img
singularity_mask bismark_genome_preparation $WORK/nextflow/work/singularity/nfcore-methylseq-1.6.1.img
singularity_mask cutadapt $WORK/singularity/trim_galore_0.6.7.sif
singularity_mask fastqc $WORK/singularity/trim_galore_0.6.7.sif
singularity_mask biopython $WORK/singularity/biopython_latest.sif python3
singularity_mask r-tidy-opt $WORK/singularity/r-tidyverse-optparse_4.1.2.sif Rscript
singularity_mask deduplicate_bismark $WORK/nextflow/work/singularity/nfcore-methylseq-1.6.1.img
singularity_mask bismark_methylation_extractor $WORK/nextflow/work/singularity/nfcore-methylseq-1.6.1.img
