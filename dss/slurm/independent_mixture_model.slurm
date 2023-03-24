#!/bin/bash

#-------------------------------------------------------
# SBATCH -J mixture_model  # sed_* will be replaced w/ sed on creation of experiment
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p development
#SBATCH -o logs/mixture_model.o
#SBATCH -e logs/mixture_model.e
#SBATCH -t 1:00:00  

#------------------------------------------------------
R_BRMS=r-brms

# override via commandline
delta=${delta:-0.25}
N=${N:-8}
brm_seed=${brm_seed:-99723}
mkdir -p model_data

# all data (maternal + paternal PO, so intragenerational data are duplicated)
output_data="model_data/dss_filtered_uncorrelated_${delta}_N${N}.bed"
# Contains just data points in the gene body
gbm_data="model_data/dss_filtered_uncorrelated_gbm_${delta}_N${N}.bed"
# Contains just data points in the gene body
nogene_data="model_data/dss_filtered_uncorrelated_nogene_${delta}_N${N}.bed"


# bed of gene body features
gbm_filter="$SCRATCH/hybrid_methylation/setup/wgbs_setup/genomes_for_bismark/Amil/cdsBoundaries.bed"


# BRMS Specific outputs
stan_fit="model_data/mixture_model_uncorrelated_fit_${delta}_N${N}.rds"
brm_draws="model_data/mixture_model_uncorrelated_draws_${delta}_N${N}.rds"
theta_quants="model_data/mixture_model_uncorrelated_theta_quantiles_${delta}_N${N}.csv"
reg_line="model_data/mixture_model_uncorrelated_regression_line_${delta}_N${N}.csv"

# BRMS Specific outputs for gbm-specific run
gbm_fit="model_data/mixture_model_uncorrelated_gbm_fit_${delta}_N${N}.rds"
gbm_brm_draws="model_data/mixture_model_uncorrelated_gbm_draws_${delta}_N${N}.rds"
gbm_theta_quants="model_data/mixture_model_uncorrelated_gbm_theta_quantiles_${delta}_N${N}.csv"
gbm_reg_line="model_data/mixture_model_uncorrelated_gbm_regression_line_${delta}_N${N}.csv"

# Non-coding run
nogene_fit="model_data/mixture_model_uncorrelated_nogene_fit_${delta}_N${N}.rds"
nogene_brm_draws="model_data/mixture_model_uncorrelated_nogene_draws_${delta}_N${N}.rds"
nogene_theta_quants="model_data/mixture_model_uncorrelated_nogene_theta_quantiles_${delta}_N${N}.csv"
nogene_reg_line="model_data/mixture_model_uncorrelated_nogene_regression_line_${delta}_N${N}.csv"



# Merge, format, and subset the data
r-tidy-opt scripts/dss_independent_merge_data.r $output_data $delta $N
# Create the GBM /noncoding filtered data
bedtools intersect -header -u -a $output_data -b $gbm_filter  > $gbm_data
bedtools intersect -header -v -a $output_data -b $gbm_filter  > $nogene_data

# run the model on all data
$R_BRMS scripts/dss_mixture_model_run.r $output_data $stan_fit $brm_seed

# run the model on GBM data
$R_BRMS scripts/dss_mixture_model_run.r $gbm_data $gbm_fit ${brm_seed}1   # Figure out how to add to it
# run the model on non-coding data
$R_BRMS scripts/dss_mixture_model_run.r $nogene_data $nogene_fit ${brm_seed}2   # Figure out how to add to it

# then, extract the outputs
$R_BRMS scripts/dss_mixture_model_extract.r $stan_fit $brm_draws $theta_quants $reg_line
$R_BRMS scripts/dss_mixture_model_extract.r $gbm_fit $gbm_brm_draws $gbm_theta_quants $gbm_reg_line
$R_BRMS scripts/dss_mixture_model_extract.r $nogene_fit $nogene_brm_draws $nogene_theta_quants $nogene_reg_line

