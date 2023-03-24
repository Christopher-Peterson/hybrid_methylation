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
delta=${delta:-0.5}
N=${N:-5}
brm_seed=${brm_seed:-16245954}

# out_dir="par_off_single_out_${delta}_N${N}"
base_data="out/joint_with_N_delta_${delta}_N${N}.rds"
# main_data="out/joint_with_N_delta_${delta}_N${N}.bed" # produced by merge_methyl_counts.slurm
# delta_data="${out_dir}/dml_with_po_delta_${delta}_N${N}.bed"


in_data="out/dss_filtered_data_delta_${delta}_N${N}.bed"
stan_fit="out/mixture_model_fit_${delta}_N${N}.rds"

# BRM Specific outputs
brm_draws="out/mixture_model_draws_${delta}_N${N}.rds"
theta_quants="out/mixture_model_theta_quantiles_${delta}_N${N}.csv"
reg_line="out/mixture_model_regression_line_${delta}_N${N}.csv"

# First, run the model
$R_BRMS scripts/dss_mixture_model_run.r $in_data $stan_fit $brm_seed

# then, extract the outputs
$R_BRMS scripts/dss_mixture_model_extract.r $stan_fit $brm_draws $theta_quants $reg_line

# Copy output files to figure data directory
fig_dat_dir=$SCRATCH/hybrid_methylation/figures/data
mkdir -p $fig_dat_dir
cp $in_data $fig_dat_dir &
cp $theta_quants $fig_dat_dir &
cp $reg_line $fig_dat_dir &

wait
