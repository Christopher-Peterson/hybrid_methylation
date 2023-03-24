# Merge all of the datasets and order them for paternal/maternal differences
suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(tidyr); library(glue)
  library(readr);
  library(rlang);
  })
if(!exists('argv')) argv = commandArgs(TRUE)

# print(argv)

delta_po_data = argv[1] %|% 'par_off_single_out/dml_with_po.bed'
out_file = argv[2] %|% 'out/dss_filtered_data.bed'
# metadata tables
parent_role_data = argv[3] %|% (Sys.getenv('SCRATCH') |> file.path('hybrid_methylation/metadata/dam_sire.csv'))
genome_order_data = argv[4] %|% 'offspring_parents.tsv'
dss_order_data = argv[5] %|% 'dss_delta_order.tsv'

### Determine Sample Metadata ####

# For the final version of this analysis, we want $\Delta$ to be Maternal - Paternal.
# For some cases, this has been lost; however, the data trail exist to identify 
# which samples need to have their signs flipped. 

parent_role = read_csv(parent_role_data, show_col_types = FALSE) |> 
  pivot_longer(DAM:SIRE, names_to = 'role', values_to = 'parent') |> 
  rename(offspring = ID) |> 
  mutate(role = recode(role, DAM='maternal', SIRE='paternal'))

genome_order = read_tsv(genome_order_data, show_col_types = FALSE,
                        col_names = c('offspring', 'g2', 'g1')) |> 
  pivot_longer(g1:g2, names_to = 'genome', names_prefix = 'g', values_to = 'parent') |> 
  mutate(offspring_asm = glue("{offspring}.{genome}"))

dss_order = read_tsv(dss_order_data, show_col_types = FALSE,
                      col_names = c('first', 'second')) |> 
  mutate(correct = glue('{first}_{second}'), 
         flip = glue("{second}_{first}")) |> 
  select(correct, flip) |> 
  pivot_longer(everything(), names_to = 'order', values_to = 'concat')

# Combine ordering info

sample_metadata = genome_order |> 
  left_join(parent_role, by = c('parent', 'offspring')) |> 
  select(-genome) |> 
  pivot_longer(c(parent, offspring_asm), names_to = 'generation', values_to = 'id') |> 
  pivot_wider(names_from = role, values_from = id) |> 
  mutate(concat = glue("{maternal}_{paternal}")) |> 
  left_join(dss_order, by = 'concat') |> 
  mutate(sign = recode(order, correct = 1L, flip = -1L))

#### Merge & create output data ####

offspring_flip = sample_metadata |> filter(order == 'flip') |> distinct(offspring) |> pull(offspring)

# read in and do some initial filtering of the data
po_data = read_tsv(delta_po_data, col_select = -c(chrom2, start2, stop2), show_col_types = FALSE) |> 
  filter(offspring==offspring2) |>  # There are false combos in this dataset; remove them
  select(-sire, -dam, -prob_5, -offspring2) |> 
  # Correct the sign
  mutate(flip_sign = if_else(offspring %in% offspring_flip, -1L, 1L),
         across(c(delta_o, delta_p), \(x) x * flip_sign)) |> 
  select(-flip_sign) |> 
  # Define delta_p - delta_o
  mutate(diff_delta = sign(delta_p) * (delta_p - delta_o), 
         diff_delta_se = sqrt(delta.se_p ^2 + delta.se_o ^ 2)) |> 
  # Add in parental role information
  left_join(parent_role, by = c('offspring', 'parent'))

# Now reorganize
out_data = po_data |> 
  relocate(1:3, offspring, parent, role, min_N,
           delta_p, delta.se_p, delta_o, delta.se_o,
           diff_delta, diff_delta_se, po_delta, po_delta_se) |> 
  arrange(`#chrom`, stop, offspring, role)
write_tsv(out_data, out_file)
# out_data |> glimpse()
# Create metadata
glue("## Metadata for ${out_file}
     #chrom, start, stop: bedgraph formatted positioning of locus. Start will always be 1 less than stop because it is zero indexed
     offspring: offspring ID
     parent: parent ID
     role: role of parent (paternal / maternal)
     min_N: lowest number of reads across paternal, maternal, and two allele-specific offspring genomes
     ## The next six lines are repeated twice, once for paternal and once for maternal
     ## Filter the data by role if they are the focus of the analysis
     delta_p: maternal - paternal methylation difference
     delta.se_p: SE of delta_p
     delta_o: Offspring maternal allele - paternal allele methylation difference
     delta.se_o: SE of delta_o
     diff_delta: Integenerational change in differential methylation; sign(delta_p) * (delta_p - delta_o)
     diff_delta_se: SE of diff_delta
     po_delta: parent - offspring methylation for maternal or paternal lineage (based on role)
     po_delta_se: SE of po_delta
     ") |> write_file(paste0(out_file, '_metadata.txt'))
  
  
### Old ####
# dam_sire = tribble(
#   ~'offspring',	~'dam',	~'sire',
#   'X1',	'A1',	'A2',
#   'X2',	'A2',	'A1',
#   'X3',	'A3',	'A4',
#   'X4',	'A4',	'A3',
#   'X5',	'A1',	'A3',
#   'X6',	'A1',	'A4',
#   'X7',	'A2',	'A3',
#   'X8',	'A2',	'A4',
# ) |> pivot_longer(dam:sire, names_to = 'role', values_to = 'parent') 
# Genome association (g1 vs g2) for allele-specific offspring loci
# This is the data in  dss/offspring_parents.tsv
# allele_specific_genomes = tribble(
#   ~'offspring', ~'g2', ~'g1', # Why did I make the order like this?
#   "X1",    "A2",      "A1",
#   "X2",    "A2",      "A1",
#   "X3",    "A4",      "A3",
#   "X4",    "A4",      "A3",
#   "X5",    "A3",      "A1",
#   "X6",    "A4",      "A1",
#   "X7",    "A2",      "A3",
#   "X8",    "A2",      "A4",
# ) |> pivot_longer(g1:g2, names_to = 'genome', names_prefix = 'g', values_to = 'parent') |> 
#   mutate(offspring_asm = glue("{offspring}.{genome}"))
# 
# dss_order = tribble(~'first', ~'second',
#                     "A1",      "A2",
#                     "A1",      "A3",
#                     "A1",      "A4",
#                     "A2",      "A3",
#                     "A2",      "A4",
#                     "A3",      "A4",
#                     "X1.1",    "X1.2",
#                     "X2.1",    "X2.2",
#                     "X3.1",    "X3.2",
#                     "X4.1",    "X4.2",
#                     "X5.1",    "X5.2",
#                     "X6.1",    "X6.2",
#                     "X7.2",    "X7.1",
#                     "X8.2",    "X8.1",
# ) |> mutate(correct = glue('{first}_{second}'), 
#             flip = glue("{second}_{first}")) |> 
#   select(correct, flip) |> 
#   pivot_longer(everything(), names_to = 'order', values_to = 'concat')
# 
# result_table = allele_specific_genomes |> 
#   left_join(dam_sire, by = c('parent', 'offspring')) |> 
#   select(-genome) |> 
#   pivot_longer(c(parent, offspring_asm), names_to = 'generation', values_to = 'id') |> 
#   pivot_wider(names_from = role, values_from = id) |> 
#   mutate(concat = glue("{dam}_{sire}")) |> 
#   left_join(dss_order, by = 'concat')
