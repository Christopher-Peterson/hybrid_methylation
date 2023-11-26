
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Analysis for “Mixed patterns of intergenerational DNA methylation inheritance in *Acropora*

*Christopher R. Peterson, Carly B. Scott, Rashin Ghaffari, Groves B.
Dixon, and Mikhail V. Matz*

## Setup

This project used Lonestar6 at the Texas Advanced Computing Center. Code
has been written specifically for this system, although it should work
on other HPC clusters with some modification.

1.  From `$SCRATCH`, clone this repository (or download the zip/tarball
    and extract it)

``` bash
cds
git clone https://github.com/Christopher-Peterson/hybrid_methylation.git
```

2.  Setup and configure the [Docker/Singularity images](setup/docker/).

Data from several stages of this pipeline are available: (details)

Contributors to this repository should also follow [these
instructions](setup/dev).

## Allele-Specific Methylation Analysis Pipeline

1.  [Get the WBGS data, trim the reads, and prepare the reference
    genome](setup/wgbs_setup/).

2.  [Align the trimmed parental reads and extract the methylation
    counts](alignment_experiments/).

3.  [Prepare and align mdRAD data](setup/md_rad/).

4.  [Variant call the parents to identify fixed SNP differences;
    generate masked genomes for allele-specific alignment; align
    offspring](snps/).

5.  [Estimate differential methylation between alleles; run Bayesian
    models to estimate inheritance](dss/).
