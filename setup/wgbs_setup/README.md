
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
tar –xvzf JA20120.tar.gz &

cd JA20120_download

# There should be 12 directories in here
ls | wc -l # 13
rm *json # Remove an unneeded json file

# move these things...
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
ls *.fastq | wc -l  # 24
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
```

## Trim the reads

### Test trimming on a data subset

First, we’re going to create a subset of the fastqs to test the trimming
on:

``` bash
NTEST=1000000 # Size of test in lines
cds hybrid_methylation/setup/wgbs_setup

mkdir -p trim_test/fastqs
# Copy the first part of each file
cd raw_data/fastqs
for file in *.fastq
do head -n $NTEST $file > ../../trim_test/fastqs/${file}
done

cd ../.. # wgbs_setup
```

Next, trim the test files.

``` bash
sbatch slurm/trim_pe_test.slurm
# This should take < 5 minutes on ls6
```

Once that’s done, run the fastqc on both sets of test data.

``` bash
mkdir trim_test/fastqc_results
sbatch slurm/fastqc_test.slurm
```

(Look at the results some time.)

### Trim the whole reads

The full fastq’s should take around 20-30 minutes to trim.

``` bash
sbatch slurm/trim_pe_full.slurm

# Check results
ls $SCRATCH/hybrid_methylation/reads/trimmed | wc -l # 24
```

They should end up in the `reads/trimmed` directory

## Prepare the Reference Genome

Most of the work for this is done by the slurm script. First, the script
files need to have their permissions enabled.

``` bash
chmod +x scripts/*py
sbatch slurm/setup_reference_genome_amilsel.sh # This one ends in .sh because it's more of a regular script that just happens to work best on slurm
```

This will create a `genomes.sh` script in the main hybrid_methylation
directory that can be sourced by other scripts.

### Create a Bismarkified version of Chromosome 1 for test mapping

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
