
# Methylation reports
# Run this from the base directory of the alignment 

suppressPackageStartupMessages({
  library(tidyverse)
  library(glue)
  library(rlang)
})
### Command Arguments
if(!exists("argv")) argv = commandArgs(TRUE)
# Arg 1: Minimum coverage for a read to be counted
min_coverage = argv[1] |> as.integer() %|% 20L
# Arg 2: Either the highest percent or count of heterogenous reads that will still be considered fixed
max_het_arg = argv[2] |> as.numeric() %|% 0.02

### Read some required metadata and file names in:
cross_data = read_csv("../../metadata/dam_sire.csv") # parent-offspring metadata
parent_files = dir("methyl_extract", pattern = "lane1-A.+cov.gz", full.names = TRUE) |> 
  set_names(c(paste0("A", 1:4)))
offspring_files = dir("methyl_extract", pattern = ".+-X.-.+cov.gz", full.names = TRUE)

### Functions ####

#' @param file bismark.cov.gz file (can also work uncompressed)
#' @param min_cov minimum coverage requirement
#' @param max_het heterogenous methylation proportions (if < 1) or counts (if >= 1) below this will still be considered fixed
#' @return a data frame of fixed parental methylation positions & calls
fixed_parental_methyl = function(file, min_cov = min_coverage, 
                                      max_het = max_het_arg) {
  name = file |> basename() |>  str_remove("^lane.-") |> str_remove("-.+$")
  # out_file = glue::glue("{out_dir}/{name}_fixedMethyl.rds")
  data = read_tsv(file = file, col_types = "ciinii", 
                  col_names = c("chrom", "start", "end", "prec_methyl", "n_meth", "n_unmeth"),
                  col_select = -end) |>
    filter(chrom != 'lambda_phage') |> 
    mutate(coverage = n_meth + n_unmeth, 
           prec_methyl = prec_methyl/100,
           minor_count = pmin(n_meth, n_unmeth),
           max_het = max_het, # Max het is either a percentage or a count; convert to a count
           max_het_count = if_else(max_het >= 1, max_het, max_het * coverage)
           ) |> 
    select(-max_het) |> 
    filter(coverage >= min_cov,
           minor_count <= max_het_count) |> 
           # prec_methyl <= fixed_thresh | prec_methyl >= (1 - fixed_thresh) ) |> 
    mutate(meth_state = as.integer(round(prec_methyl))) |> 
    select(chrom, start, meth_state, coverage, minor_count) |>
    mutate(ID = name)
  # write_rds(data, out_file)
  # invisible()
  data
}

#' Combine two parental fixed methyl points (only keep overlaps)
#' @param parent_data named list of fixed parental methylation sites
#' @param dam ID of mother
#' @param sire ID of father
#' @return combined data frame with all positions that both parents are fixed for methylation
joint_methyl_parents = function(parent_data, dam, sire) {
  dam_data = parent_data[[dam]] |> 
    #read_rds(glue("methyl_extract/fixed_loci/{dam_ID}_fixedMethyl.rds")) |> 
    rename(dam_meth = meth_state, dam = ID, dam_coverage = coverage, dam_minor = minor_count)
  sire_data = parent_data[[sire]] |> 
    # read_rds(glue("methyl_extract/fixed_loci/{sire_ID}_fixedMethyl.rds")) |> 
    rename(sire_meth = meth_state, sire = ID, sire_coverage = coverage, sire_minor = minor_count)
  # browser()
  # columns: chrom, start, meth_state, ID
  inner_join(dam_data, sire_data, by = c("chrom", "start")) #|> 
    # filter(dam_meth != sire_meth)
  # contrasting_methyl
}

#' Read offspring methlation calls and filter to include only fixed parental positions
#' @param file path to a bismark.cov.gz file
#' @param parent_data named list of all fixed parental positions 
#' @param corssing data frame listing the DAM and SIRE for each offspring ID
#' @param min_cov minimum coverage for offspring loci
#' @return a data frame of offspring methylation percent at fixed parental loci
filter_offspring_methyl = function(file, parent_data, crossing = cross_data,
                                   min_cov = min_coverage#, out_dir = "methyl_extract/fixed_loci"
                                   ) {
  # browser()
  name = file |> basename() |>  str_remove("^lane.-") |> str_remove("-.+$")
  # browser()
  parents = crossing |> filter(ID == name)
  joint_parents = joint_methyl_parents(parent_data, parents$DAM, parents$SIRE)
  # out_file = glue::glue("{out_dir}/{name}_contrasts.rds")
  offspring_data = read_tsv(file = file, col_types = "ciinii", col_names = c("chrom", "start", "end", "prec_methyl", "n_meth", "n_unmeth"), col_select = -end) |>
    filter(chrom != 'lambda_phage') |> 
    mutate(coverage = n_meth + n_unmeth, ID = name) |> 
    filter(coverage >= min_cov)
  result = joint_parents |> left_join(offspring_data, by = c("chrom", "start")) |> 
    filter(!is.na(prec_methyl))
  # write_rds(result, out_file)
  # invisible()
  result
}

### Run the functions ####

parent_data_list = map(parent_files, fixed_parental_methyl)

parent_species = tribble(
  ~"ID", ~"species",
  "A1", "Asel",
  "A2", "Asel",
  "A3", "Amil",
  "A4", "Amil",
)

offspring_data = map_dfr(offspring_files, filter_offspring_methyl, parent_data = parent_data_list) |> 
  left_join(parent_species |> rename(dam = ID, dam_spp = species)) |> 
  left_join(parent_species |> rename(sire = ID, sire_spp = species)) |> 
  mutate(species = if_else(dam_spp == sire_spp, dam_spp, "hybrid"),
         prec_methyl = prec_methyl / 100,
         parent_status = case_when(
           dam_meth == sire_meth & dam_meth == 1L ~ "Both Parents Methylated",
           dam_meth == sire_meth & dam_meth == 0L ~ "Both Parents Unmethylated",
           dam_meth == 1L ~ "Mother Methylated, Father Unmethylated", 
           dam_meth == 0L ~ "Mother Unmethlated, Father Methylated")
  )

## Save the Output ####
out_dir = "methyl_extract/fixed_loci"
dir.create(out_dir, showWarnings = FALSE)
out_file = glue("{out_dir}/results_{min_coverage}_{max_het_arg}.rds")

list(result = offspring_data, args = c(min_coverage, max_het_arg)) |> write_rds(out_file)

#  # This code jumps into the singularity R
# module load tacc-singularity
# HDIR=/home
# IMAGE=/work/04386/crpeters/ls6/singularity/r-tidyverse-optparse_4.1.2.sif
# CMD=Rscript
# singularity shell -H $HDIR $IMAGE $CMD "$@"
# R
