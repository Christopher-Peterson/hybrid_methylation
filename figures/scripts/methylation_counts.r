# count_beds/joint
library(tidyverse)
library(glue)
library(ggtext)
read_bsseq_cov = \(file) {
  dat = read_tsv(file, col_names = c('chr', 'pos','end', 'meth_prec', 'N_meth', 'N_unmeth'), 
                 col_select = -c(end, meth_prec))
  # browser()
  dat |> dplyr::mutate(N = N_meth + N_unmeth) |> 
    dplyr::select(chr, pos, N, X = N_meth)
}


adult_list = local({
  adult_names = c("A1", "A2", "A3", "A4")
  dir("parents", pattern = "*A[1-4].*cov.gz", full.names = TRUE) |> 
    set_names(adult_names)
})
offspring_list = local({
  files = dir("offspring", pattern = "X[1-8].genome[1|2].bismark.cov.gz", full.names = TRUE) 
  nms = files |> basename() |> str_remove('.bismark.cov.gz') |> str_remove("genome")
  files |> set_names(nms)
})

offspring_parents = read_tsv('offspring_parents.tsv', col_names = c('cross', 'par2', 'par1')) |> 
  mutate(off2 = paste0(cross, '.2'), off1 = paste0(cross, '.1'))

dir.create('intersected_counts')
# These two functions make 
intersect_bsseq_cov = \(cross, par1, par2, off1, off2, min_reads = 5, max_reads = 400) {
  p1 = read_bsseq_cov(adult_list[par1])
  o1 = read_bsseq_cov(offspring_list[off1])
  join1 = inner_join(p1, o1, by = c('chr', 'pos'), suffix = c('_p1', '_o1'))
  rm(p1, o1)
  p2 = read_bsseq_cov(adult_list[par2])
  o2 = read_bsseq_cov(offspring_list[off2])
  join2 = inner_join(p2, o2, by = c('chr', 'pos'), suffix = c('_p2', '_o2'))
  rm(p2, o2)
  # browser()
  full_join = 
    inner_join(join1, join2, by = c('chr', 'pos'))
  count_file = bind_rows(join1 |> mutate(genome = 'genome1'),
            join2 |> mutate(genome = 'genome2'),
            full_join |> mutate(genome = 'full')) |>
    mutate(N = pmin(N_p1, N_p2, N_o1, N_o2, na.rm = TRUE)) |> 
    group_by(genome, chr, N) |> 
    summarize(count = n())
  filtered_join = full_join |> 
    filter(if_all(c(N_p1, N_p2, N_o1, N_o2), 
                   list(\(x) x >= min_reads, \(x) x <= max_reads)))
  out_file = file.path('intersected_counts', paste0('cross_', cross, '.tsv'))
  filtered_join |> rename(`#chrom` = chr) |> 
    mutate(cross = cross) |> write_tsv(out_file)
  count_file |> mutate(cross = cross)  
}
methyl_count_report = function(file_pair) {
  p1 = read_bsseq_cov(file_pair[1])
  p2 = read_bsseq_cov(file_pair[2])
  # browser()
  joint_data = 
    inner_join(p1, p2, by = c('chr', 'pos'), suffix = c('_1', '_2')) |> 
    mutate(N_min = pmin(N_1, N_2))
  nms = names(file_pair)
  count_data = joint_data |> mutate(N = N_min, name = paste(nms, collapse = " x ")) |> 
    bind_rows(p1 |> mutate(name = nms[1] ),
              p2 |> mutate(name = nms[2])
    ) |> group_by(name, N) |> 
    summarize(count = n())
  # browser()
  count_data
}

library(parallel)
offspring_pairs = lapply(1:8, \(x) glue("X{x}.{c(1:2)}"))
offspring_summary = map(offspring_pairs, \(x) offspring_list[x]) |> 
  map(methyl_count_report) |> bind_rows()
adult_summary = combn(paste0("A", 1:4), 2, simplify = FALSE) |> 
  map(\(x) adult_list[x]) |>
  map(methyl_count_report) |> bind_rows()

full_cross_counts = offspring_parents |> 
  pmap_dfr(intersect_bsseq_cov) |> 
  group_by(cross, genome) |>
  select(-chr) |> # summarize(count = sum(count)) |>
  uncount(count) |> 
  chop(N) |> 
  mutate(dens = lapply(N, \(x) {
    # browser();
    dens = density(log10(x), from = 0)
    tibble(dens$x, dens$y)
  })) |> unpack(dens)

full_cross_counts |> write_rds('out/full_cross_methyl_counts.rds', compress = 'bz2')
# full_cross_counts  = read_rds('out/full_cross_methyl_counts.rds')

pairwise_methyl_counts =  bind_rows(offspring_summary, adult_summary)
  pairwise_methyl_counts |> write_rds('out/pairwise_methyl_counts.rds')
  
  
### Switch to local   ####
pairwise_methyl_counts =  read_rds('figures/data/pairwise_methyl_counts.rds')
# full_cross_counts = read_rds('figures/data/full_cross_methyl_counts.rds')
full_cross_counts_smry = full_cross_counts |>
  select(-dens) |> 
  unchop(N) |> 
  group_by(genome, cross, N) |>
  summarize(count = n()) |>
  ungroup() |> 
  chop(c(N, count)) |> 
  mutate(group = recode(genome, full = 'Complete', genome1 = 'Lineage 1', genome2 = 'Lineage 2'),
         overlap = if_else(genome == 'full', 'Complete', 'Specific Allele (Parent / Offspring)')) |> 
  select(cross, overlap, group, N, count) 
pairwise_methyl_counts_details = pairwise_methyl_counts |> 
  mutate(class = substr(name, 1, 1) |> recode(X = 'Offspring', A = 'Parent'),
         paired = str_detect(name, ' x ') |> if_else('Pair', 'Single')) |> 
  chop(c(N, count))

individual_level_counts = offspring_parents |> 
  pivot_longer(-cross, names_to = 'ind', values_to = 'name') |> 
  left_join(pairwise_methyl_counts_details) |> 
  mutate(overlap = recode(ind, par1 = "One Parent", par2 = "One Parent", 
                          off1 = 'One Offspring Allele', off2 = 'One Offspring Allele')) |> 
  select(cross, overlap, group = ind, N, count)
within_generation_pairs = offspring_parents |> 
  mutate(parents = paste(pmin(par1, par2), pmax(par1, par2), sep = ' x '),
         offspring = paste(off1, off2, sep = ' x ')) |> 
  select(cross, parents, offspring) |> 
  pivot_longer(-cross, names_to = 'pair', values_to = 'name') |> 
  left_join(pairwise_methyl_counts_details) |> 
  mutate(overlap = recode(pair, parents = 'Both Parents', offspring = 'Both Offspring Alleles'),
         group = overlap) |> 
  select(cross, overlap, N, count, group)

pretty_parents = \(x) recode(x, 
                             A1 = 's<sub>1</sub>', A2 = 's<sub>2</sub>',
                             A3 = 'm<sub>1</sub>', A4 = 'm<sub>2</sub>')

offspring_details = tribble(~"offspring",  ~"maternal", ~"paternal",
                            "X1",      "A1",      "A2",
                            "X3",      "A3",      "A4",
                            "X5",      "A1",      "A3",
                            "X6",      "A1",      "A4",
                            "X7",      "A2",      "A3",
                            "X8",      "A2",      "A4",
                            "X2",      "A2",      "A1",
                            "X4",      "A4",      "A3") |> 
  mutate(across(maternal:paternal, pretty_parents),
         pretty_offspring = glue("{maternal} × {paternal}")
  ) |> select(cross = offspring, pretty_offspring) |> 
  arrange(cross) |> 
  mutate(pretty_offspring = fct_inorder(pretty_offspring))

complete_counts = bind_rows(full_cross_counts_smry, individual_level_counts, within_generation_pairs)  |> 
  mutate(size = map_int(N, length)) |> 
  arrange(desc(cross), desc(size))  |> 
  mutate(overlap = if_else(overlap == 'Specific Allele', 'Specific Allele (Parent / Offspring)', overlap)) |> 
  mutate(overlap = fct_inorder(overlap)) |> 
  unchop(c(N, count)) |> 
  left_join(offspring_details, 'cross') 

# complete_counts |> 
#   ggplot(aes(x = N, y = count, color = overlap, group = group)) + 
#   geom_line() + 
#   scale_y_log10('Number of CpG Sites') + scale_x_log10('Number of Reads') + 
#   scale_color_viridis_d(option = 'turbo') + 
#   theme_classic() + 
#   facet_wrap(~pretty_offspring, nrow = 2) + 
#   theme(strip.background = element_blank(), strip.text = element_markdown())
#   

methylation_counts_figure = complete_counts |> 
  filter(N <= 100) |> 
  ggplot(aes(x = N, y = count, color = pretty_offspring, group = paste0(group, pretty_offspring))) + 
  geom_line() + 
  scale_y_continuous('Number of CpG Sites', trans = 'log10', breaks = 10^(0:8),
                labels = scales::label_number(scale_cut = scales::cut_si(''))) + 
  scale_x_log10('Minimum number of reads at CpG Site') + 
  scale_color_viridis_d('Cross', option = 'turbo', end = .9) + 
  theme_classic() + 
  facet_wrap(~overlap, scales = 'free_x') + 
  theme(strip.background = element_blank(), strip.text = element_markdown(),
        legend.text = element_markdown())
ggsave('figures/methylation_counts_by_overlap.png',methylation_counts_figure, dpi = 300, width = 12, height = 8)
  
