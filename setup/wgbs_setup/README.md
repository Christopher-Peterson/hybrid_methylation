
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Whole Genome Bisulfite Data setup

This section documents acquiring the WBGS reads, trimming them, and
preparing the reference genome.

## Acquire the data files

Copy the files from their archive into a raw data folder and extract it.

``` bash
cds hybrid_methylation/setup/wgbs_setup
# Setup raw_data and slurm sub-directories
mkdir logs jobs raw_data raw_data/fastqs
cd raw_data

# Copy data from Grove's archive
cp  /corral-repl/utexas/tagmap/dixon_backups/file_sharing_test_dir/JA20120.tar.gz .

# Untar
tar –xzf JA20120.tar.gz 

cd JA20120_download

# There should be 12 directories in here
ls | wc -l # 13
rm *json # Remove an unneeded json file
```

Each of these directories contains two `fasta.gz` files (the paired
ends); move them to a new, single directory and unzip them.

``` bash
mv */*gz ../fastqs
cd .. # raw_data
ls fastqs/*gz | wc -l # should be 24
# Unzip in parallel on a dev node
cd .. # wgbs_setup
sbatch slurm/unzip_fastq.slurm
```

When that’s done, check that the counts are correct.

``` bash
# Expected numbers in comments:
ls raw_data/fastqs/*.fastq | wc -l  # 24
ls lane1*.fastq | wc -l # 12
ls lane2*.fastq | wc -l # 12
ls *_R1_*.fastq | wc -l # 12
ls *_R2_*.fastq | wc -l # 12
```

Create a symbolic link in the main `$SCRATCH` directory for easier
access.

``` bash
cds hybrid_methylation
mkdir reads
cd reads
ln -s ../setup/wgbs_setup/raw_data/fastqs raw

cds hybrid_methylation/setup/wgbs_setup
```

## Trim the reads

### Trimming with Trim Galore

We’re using Trim Galore.

``` bash
sbatch slurm/trim_galore_pe_full.slurm

# This produces trim files with the wrong name; adjust them to have the same names as the non trim-galore method previously tried 
# Need to go from lane1-A1-A_S1_L001_R1_001_val_1.fq  
  # to lane1-A1-A_S1_L001_1.trim

function rename_file {
  local a=$1
  local b=${a/_R?_001_val/}
  local c=${b/.fq/.trim}
  mv $a $c
}
for fq in raw_data/trimmed/*fq; do
    rename_file $fq ; done
```

Finally, we want to create a link to the trimmed files in the main
directory

``` bash
trim_dir=$PWD/raw_data/trimmed
link_dir=$SCRATCH/hybrid_methylation/reads/trimmed
mkdir -p $link_dir
for fq in $trim_dir/lane*trim; do
  ln -s $fq $link_dir/$(basename $fq)
done
```

## Prepare the Reference Genome

This acquires the Amil v. 2.1 reference genome and its GFF (may need to
be altered to target an external source), runs the bismark genome script
on it, and generates feature windows. Most of the work for this is done
by the slurm script. First, the script files need to have their
permissions enabled.

``` bash
chmod +x scripts/*py
sbatch slurm/setup_reference_genome_amil.sh # This one ends in .sh because it's more of a regular script that just happens to work best on slurm
```

This will create a `genomes.sh` script in the main hybrid_methylation
directory that can be sourced by other scripts.

### Create a Bismarkified version of Chromosome 1 for test mapping

*Not required for replication.*

For testing purposes, we will also create a version of the reference
genome that’s only for chromosome 1. First, we need to split the genome
into individual chromosomes. Version 2.1 of the *A. mil* genome has 14
chromosomes and 840 unplaced elements. We will only focus on the
chromosomes.

``` bash
####### new Script
# Create the directory for the chromosomes
mkdir $SCRATCH/hybrid_methylation/genomes/Amil_v2.1_chrs
cd $SCRATCH/hybrid_methylation/genomes/Amil_v2.1_chrs

GENOME_FASTA=../Amil_v2.1/GCF_013753865.1_Amil_v2.1_genomic.fna 

# Use csplit to split the genome file into 
# context-determined pieces, with new filenames xx00, xx01, ...
# -s  silent
# -z Suppress generation fo zero-length outputs
# This treats both />/ and {*} as splitting patterns

csplit -s -z $GENOME_FASTA  '/>/' '{*}'

# The first 14 files are the chromosomes
for i in {xx0?,xx1[0-3]} # first 14 chromosomes
do 
  # Extract the chromosome number
  n=$(head $i -n1 | sed 's/.*chromosome //; s/,.*//')
  # Rename chromosomes
  mv "$i" "chr_$n.fasta"
done
# Move the unplaced scaffold files
mkdir unplaced_scaffolds
mv xx* unplaced_scaffolds
# Renamed these to their scaffold numbers
cd unplaced_scaffolds
for i in xx* 
do 
  name=$(head $i -n1 | sed 's/.*Amil_v2.1 //; s/,.*//')
  mv "$i" "$name.fasta"
done
cd ..
```

Now, we’ll run the bismark prep script on just chromosome 1.

``` bash
# Copy chr1 to new directory
mkdir ../Amil_chr1
cp chr_1.fasta ../Amil_chr1

# Convert to bismark style genome
SLURM_DIR=$SCRATCH/hybrid_methylation/setup/wgbs_setup/slurm
sbatch $SLURM_DIR/setup_split_chrom_genome.sh
```
