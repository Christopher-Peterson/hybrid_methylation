
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Tag-seq Data setup

This section documents acquiring the tag-seq reads, trimming, mapping,
and counting them.

## Acquire the data files

Copy the files from their archive into a raw data folder and extract it.

``` bash
cds hybrid_methylation/tagseq_setup
# Setup raw_data and slurm sub-directories
mkdir logs jobs raw_data fastqs

# Copy data from Carly's work (TODO: MAKE SURE THESE ARE ARCHIVED SOMEWHERE)
cp  /work/06909/cbscott/JA21076_download/fastqs/*tseq*.fastq.gz ./fastqs

ls fastqs/*gz | wc -l # should be 72
# Unzip in parallel on a dev node
sbatch slurm/unzip_fastq.slurm
```

When that’s done, check that the counts are correct.

``` bash
# Expected numbers in comments:
ls fastqs/*.fastq | wc -l  # 72
ls fastqs/*L001*.fastq | wc -l # 36 lane 1
ls fastqs/*L002*.fastq | wc -l # 36 lane 2
```

