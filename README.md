
<!-- README.md is generated from README.Rmd. Please edit that file -->

# The Hybrid Methylation Project

(At some point, I’ll put a real introduction here)

## Setup

This project uses Lonestar6 on TACC; I will be assuming you have access
to the system.

1.  From `$SCRATCH`, clone this repository. You may need to setup a
    Github Personal Access Token to do that on TACC.

``` bash
cds
git clone https://github.com/Christopher-Peterson/hybrid_methylation.git
```

2.  Setup and configure the [Docker/Singularity images](setup/docker/).

3.  Setup [Integrative Genomics Viewer (IVG)](setup/igv) and associated
    tools for use on TACC.

Contributors to this repository should also follow [these
instructions](setup/dev)..

## Whole Genome Bisulfide Sequence Data

1.  [Get the Data, trim the reads, and prepare the reference
    genome](setup/wgbs_setup/).

2.  [Align the trimmed reads and extract the methylation
    information](alignment_experiments/).

3.  (I am working on this part)
