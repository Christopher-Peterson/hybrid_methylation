
# Differential Methylation Analysis

## Calculate differnetial methylation

This is based off of the Bioconductor package DSS, which uses smoothing
windows & shrinkage to estimate beta-binomial dispersion parameters for
differential methylation estimates.

First, set everything up:

``` bash
cds hyb*/snps/with_md
ln -s ../../dss dss
cd dss
mkdir -p logs jobs
cp ../offspring_parents.tsv .
```

### Run DSS

This script runs per-locus DSS on all loci detected in both alleles:

``` bash
cds hyb*/snps/tg_md/dss
sbatch slurm/independent_dss.slurm
```

## Bayesian Analyses

These models are described in the publication and supporting information

The mixture model needs to be run first, as it also does setup that is
important for the other two models:

``` bash
sbatch slurm/independent_mixture_model.slurm
```

Run the bivariate model:

``` bash
sbatch slurm/run_bivar_model.slurm
```

Run the locus-specific model:

``` bash
sbatch slurm/run_consistency_model.slurm
```
