#!/bin/bash
#-------------------------------------------------------
# SBATCH -J consistency_model  # sed_* will be replaced w/ sed on creation of experiment
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p development
#SBATCH -o logs/consistency_model.o
#SBATCH -e logs/consistency_model.e
#SBATCH -t 1:00:00  

#------------------------------------------------------
R_BRMS=r-brms

# override via commandline
delta=${delta:-0.25}
N=${N:-8}
brm_seed=${brm_seed:-18245}

# out_dir="par_off_single_out_${delta}_N${N}"
# mkdir -r model_data


# Created by the dss_independent_merge_data.r script
base_data="model_data/dss_filtered_uncorrelated_${delta}_N${N}.bed"
gbm_data="model_data/gbm_loci_uncorrelated_${delta}_N${N}.bed"
stan_fit="model_data/consistency_fit_${delta}_N${N}.rds"
stan_smry="model_data/consistency_fit_${delta}_N${N}.rds"

# run the model
$R_BRMS scripts/run_consistency_model.r $base_data $gbm_data $stan_fit $stan_smry $brm_seed
    
# # Copy output files to figure data directory
# fig_dat_dir=$SCRATCH/hybrid_methylation/figures/data
# mkdir -p $fig_dat_dir
# cp $in_data $fig_dat_dir &
# cp $theta_quants $fig_dat_dir &
# cp $reg_line $fig_dat_dir &
# 
# wait
# # 