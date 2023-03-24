
<!-- README.md is generated from README.Rmd. Please edit that file -->

# md-RAD Data setup

This section documents acquiring the md-RAD reads, trimming them, and …?

## Acquire the data files

Copy the files from their archive into a raw data folder

``` bash
mkdir -p $SCRATCH/hybrid_methylation/setup/md_rad/raw_data
cds hybrid_methylation/setup/md_rad/raw_data

remote_dir=/work/06909/cbscott/JA21076_download/fastqs

# We only want the md-RAD data that is from Adults (A-*md?_S*fastq.gz) and single larvae (SL_*fastq.gz)

cp $remote_dir/A-*md?_S*fastq.gz . &

cp $remote_dir/SL-*md?_S*fastq.gz . &
```

## Prep and trim the reads for Adults

Next, unzip and concatenate the paired files. To begin with, we’ll just
do the adults

``` bash
idev
# idev-skx
for f in *gz; do pigz -d $f &
done

cd ..
mkdir concat
 
function concat_paired_reads {
  local left=$1
  local right=${left/_L001_/_L002_}
  local bn=`basename $left`
  local out=concat/${bn/_L001_R1_001/}
  cat $left $right > $out
}

for fl in raw_data/*L001_R1_001.fastq; do concat_paired_reads $fl & 
done
```

Run the preliminary filtering script on each concatenated read set

``` bash
# Expected numbers in comments:
mkdir -p scripts slurm jobs logs
chmod +x scripts/*
mkdir -p filt0 pe_filt0
idev

FILT0_LOG=logs/filt0.o

function run_filter0 { # SE version
  local IN_FILE=$1
  local OUT_FILE=${IN_FILE/.fastq/_filt0.fastq}
  # This script doesn't like changing output directories, so we move the file after the fact
  biopython scripts/filter_methylGBS_reads.py -r1 $IN_FILE -o1 $OUT_FILE >> $FILT0_LOG
  
  mv $OUT_FILE filt0
}
# PE version; this one doesn't work
function run_filter0_pe { 
  local R1=$1
  local O1=${R1/raw_data/pe_filt0}
  local R2=${R1/L001/L002}
  local O2=${R2/raw_data/pe_filt0}
  biopython scripts/filter_methylGBS_reads.py -r1 $R1 -r2 $R2 -o1 $O1 -o2 $O2 >> $FILT0_LOG
}

> $FILT0_LOG
#for file in raw_data/*S43*L001*fastq; do # PE version; doesn't work
for file in concat/S*; do  # SE version
  run_filter0 $file &
#  run_filter0_pe $file &
done

#> $FILT0_LOG
#for file in raw_data/*S43*L001*fastq; do  run_filter0_pe $file; done; cat $FILT0_LOG


# exit
# Now run the proper adapter trimming
sbatch slurm/md_cutadapt
```

## Align reads?

First, setup & Bowtie index the genome

``` bash
# Copy over the amil genome
mkdir genome
cd genome

# If it already exists on Scratch
# ln -sf ../../../genomes/Amil_v2.1/GCF_013753865.1_Amil_v2.1_genomic.fna Amil.fasta

# if Not:
# cp $STOCKYARD/tagmap-share/genomes/Amil_v2.1/GCF_013753865.1_Amil_v2.1_genomic.fna Amil.fasta

# Run the bowtie2-build indexer
idev
bowtie2-build -f --threads 20 Amil.fasta Amil
exit
cd ..
```

Run the alignment & postprocess it.

``` bash
sbatch slurm/md_align.slurm

# post-process it
sbatch slurm/md_align_post.slurm
```

Mapping for these is only about 30% unique (not ideal).
