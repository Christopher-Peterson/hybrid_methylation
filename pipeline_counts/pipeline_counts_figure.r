# pipeline counts
library(tidyverse)
library(ggtext)

files =   dir('pipeline_counts', full.names = TRUE, pattern = '.tsv') 
raw_pc = files |> 
  set_names(basename(files) |> str_remove('pipeline_') |> str_remove('counts_') |> str_remove('.tsv')) |> 
  map(read_tsv)
bind_rows(raw_pc)

raw_pc |> names()

dedup_dat = raw_pc[c('dedup', 'offspring_dedup')] |> bind_rows()
bam_dat = raw_pc[c('bam', 'offspring_bam')] |> bind_rows()
bam_dat |> distinct(column)
dedup_dat |> distinct(column)
raw_pc$trim |> distinct(column)

input_opts = c("Sequence pairs analysed in total:", "Total number of alignments analysed")
output_opts = c("Number of paired-end alignments with a unique best hit:","Total count of deduplicated leftover sequences")

conversion_table = tribble(~'step', ~'column',
  'Raw Reads', "Sequence pairs analysed in total:",
  'Aligned', "Total number of alignments analysed",
  'Deduplicated', "Total count of deduplicated leftover sequences",
  'Split (genome 1)', 'genome_1',
  'Split (genome 2)', 'genome_2'
) |> mutate(step = fct_inorder(step))

parent_ids = \(x) recode(x, A1 = "s<sub>1</sub>", A2 = "s<sub>2</sub>", A3 = "m<sub>1</sub>", A4 = "m<sub>2</sub>")

id_key = read_csv('metadata/dam_sire.csv') |> 
  arrange(ID) |> 
  mutate(Name = paste(parent_ids(DAM), parent_ids(SIRE), sep = " × "),
         generation = 'Offspring') |> 
  select(ID, Name, generation) |> 
  bind_rows(tibble(ID = paste0("A", 1:4)) |> mutate(Name = parent_ids(ID), generation = 'Parents')) |> 
  mutate(Name = fct_inorder(Name))


pipeline_counts_df = raw_pc[c('dedup', 'offspring_dedup', 'bam', 'offspring_bam', 'offspring_snpsplit')] |> bind_rows() |> 
  right_join(conversion_table, 'column') |> 
  select(-column) |> 
  rename(ID = `#ID`) |> left_join(id_key, by = 'ID')

theme_set(theme_classic())

pipeline_counts_figure = pipeline_counts_df |> ggplot(aes(x = Name, y = count)) + 
  geom_line(aes(color = step, group = step)) + 
  scale_color_viridis_d(end = .85, option = 'inferno') + 
  facet_grid(.~generation, scale = 'free_x', space = 'free_x',switch = 'x') + 
  scale_y_continuous("Read Counts", trans = 'log10', n.breaks = 12,
                     labels = scales::label_number(scale_cut = scales::cut_si(''))) +
  xlab("") + 
  theme(axis.text.x = element_markdown(),
        legend.position = c(.85, .3),
        strip.placement = 'outside', strip.background = element_blank())
ggsave('figures/pipeline_counts.png', pipeline_counts_figure, dpi = 300, width = 8, height = 5)
