# dss pairwise_post
# 
# Alternate idea:
# Convert output to a bed file
# and use bedtools for all the subsetting
# One bedfile for each FDR level
# 

# 
library(tidyverse)
library(DSS)
library(glue)
library(IRanges)
library(GenomicRanges)


offspring_info = read_tsv('offspring_parents.tsv', col_names = c('off', 'p1', 'p2')) |> 
  mutate(pair = 1:n()) |> 
  rowwise() |> 
  mutate(type = case_when(
    all(c_across(p1:p2) %in% c("A1", "A2")) ~ "A. selago F1",
    all(c_across(p1:p2) %in% c("A3", "A4")) ~ "A. milepora F1",
    TRUE ~ "Hybrid F1"
  )) |> ungroup()



read_as_df = function(files, type, pair){
  map_dfr(files, \(x) read_rds(x) |> unclass() |> structure(class = 'data.frame')) |> 
    mutate(type = type, pair = pair)
}
read_parent_data = \(p1, p2, type, pair, ...) {
  dir(base_dir, pattern = glue("({p1}|{p2})_({p2}|{p1})_chr.*.rds"), full.names = TRUE) |> 
    read_as_df(type, pair)
}
read_offspring_data = \(off, type, pair, ...) {
    base_dir |> dir(off, full = TRUE) |>  read_as_df(type, pair)
}

all_parent_dml = pmap_dfr(offspring_info, read_parent_data) 
all_offspring_dml = pmap_dfr(offspring_info, read_offspring_data) 

all_parent_dml |> write_rds("pairwise_out/dml_test/all_parents.rds")
all_offspring_dml |> write_rds("pairwise_out/dml_test/all_offspring.rds")

as_dml = function(x) structure(x, class = c("data.frame", "DMLtest"))

# bed_dir = "/work/04386/crpeters/tagmap-share/genomes/Amil_v2.1_ForBismark"

# bed_files = dir(bed_dir, ".bed", full.names = TRUE)

# Let's just run all of them...

# SO the idea here:
# Read in a bed file
# filter the parent data
# finish the script
# dir.create('window_out', recursive = TRUE)

parent = read_rds("pairwise_out/dml_test/all_parents.rds") #|> 
  # select(chr, pos, type, pair, diff, diff.se, fdr) |> as_tibble()

# parents_gr = read_rds("pairwise_out/dml_test/all_parents.rds") |> 
#   with(GRanges(chr, IPos(pos), "*", tibble(diff, diff.se)))

to_bed = \(x) x |>   select(chr, start = pos, end = pos, type, pair, chr, diff, diff.se)


all_offspring_dml = read_rds("pairwise_out/dml_test/all_offspring.rds") |> to_bed()
write_tsv(all_offspring_dml, 'bed_filter/offspring.bed', col_names = FALSE)

# New goal: 
# filter the parents to various levels and write out to bed

make_parent_bed = \(dmr_delta, fdr_level = 0.05) {
  out_name = glue::glue("bed_filter/parents_{dmr_delta}_{fdr_level}.bed")
  bed_file = parent |> 
    as_dml() |>
    callDML(delta = dmr_delta) |> 
    as_tibble() |>
    filter(fdr <= fdr_level) |> 
    to_bed()
    write_tsv(bed_file, out_name, col_names = FALSE)
  }

dmr_delta_opts = c(.2, .3, .4, .5, .55, .6, .65, .7)
map(dmr_delta_opts, make_parent_bed)


# window_analysis = \(bed_file, dmr_delta = .6, fdr_set = .05) {
#   # browser()
#   raw_input = read_tsv(bed_file, col_names =
#                          c("chr", "start", "end", "annotation")) |> 
#     filter(!is.na(start), !is.na(end))
#   ir = with(raw_input, IRanges(start = start + 1L, end = end))
#   gr = GRanges(raw_input$chr, ir, "*", mcol = raw_input |>
#                  dplyr::select(annotation) )
#   
#   overlap = findOverlaps(gr, parents_gr)
#   subset_parents = parents_gr[subjectHits(overlap)] |> as_tibble()
#   rm(ir, gr, overlap)
#   filtered_dml = subset_parents |> as_dml() |>
#     callDML(delta = dmr_delta) |> 
#     as_tibble() |>
#     filter(fdr <= fdr_set)
#   rm(subset_parents)
#   joint_dml = filtered_dml |> 
#     as_tibble() |>
#     select(type, pair, chr=seqnames, pos=start, diff, diff.se) |> 
#     inner_join(all_offspring_dml ,
#                by = c('type', 'chr', 'pair', 'pos'), suffix = c('_p', '_o'))
#   
#   
#   out_name = glue::glue("window_out/{z}_{dmr_delta}_{fdr_set}.rds",
#                         z = bed_file |> basename() |> str_remove('.bed') )
#   print(glue("writing {out_name}..."))
#   write_rds(joint_dml, out_name)
# }
# 
# bed_files |> walk(window_analysis, dmr_delta = .6)
# bed_files |> walk(window_analysis, dmr_delta = .5)
# bed_files |> walk(window_analysis, dmr_delta = .55)
# 
# 
# 
# 
# 
# dmr_delta = 0.65
# fdr_set = 0.05
# 
# # If this has trouble, re-run on full node
# parent_dml_call_filt = all_parent_dml |> as_dml() |> callDML(delta = dmr_delta) |> 
#   as_tibble() |> filter(fdr <= fdr_set)
# 
# 
# 
# joint_dml = parent_dml_call_filt |> as_tibble() |> select(type, pair, chr, pos, diff, diff.se) |> 
#   inner_join(all_offspring_dml |> select(type, pair, chr, pos, diff, diff.se),
#              by = c('type', 'chr', 'pair', 'pos'), suffix = c('_p', '_o'))
# 
# write_rds(joint_dml, glue('joint_full_{dmr_delta}.rds'))
# 
# 
# 
# 
# ### TODO:
# ### Check the directions on the non-hybrid F1's, given some of the upper-left points
# ### Make rough outline of intro, methods
# ### Subset to only include gene windows...
# ### GO analysis?  Look at gene windows?
# 
# 
# 
# parent_dmr_call  = parent_dml_call_filt |> 
#   unite('chr', sep = '/', type, pair, chr) |> 
#   as_dml() |> 
#   callDMR(delta = dmr_delta) |> 
#   as_tibble() |> 
#   separate(chr, into = c('type', 'pair', 'chr'), sep = '/')
# write_rds(parent_dmr_call, 'joint_parent_dmr.rds')
# 
# str#### Switch to local ####
# joint_dml = read_rds('dss/joint_full_0.6.rds')
# 
# 
# joint_dml |> ggplot(aes(diff_p, diff_o, color = as.factor(pair |> pmin(5)))) + 
#   facet_wrap(~type) + 
#   geom_abline(slope = 1, intercept = 0, linetype = 2) + 
#   geom_abline(slope = 0, intercept = 0, linetype = 1) + 
#   geom_point() + coord_fixed() + theme_classic() + 
#   scale_color_brewer(type = 'qual', palette = "Set1") +
#   # scale_color_viridis_c("Offspring\nDifference SE") +
#   xlab("Parental Methylation Difference") + 
#   ylab("F1 Methylation Difference")
# 
# # So the next idea:
# 
# # Filter the offspring dml into ones that are clearly > 0.5 and those that are < 0.2
# # For many of them, there seems to be 'chunks'
# # Then look at GO terms?  
# # 