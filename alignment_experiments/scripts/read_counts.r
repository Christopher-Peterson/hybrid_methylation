# Read a .wig file (of genome counts, created by IGV)
# Convert it to a useful data frame

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(rlang) # Necessary for processing command args
})

# Read Command Args ####
if(!exists("argv")) argv = commandArgs(TRUE)
in_file  = argv[1]   %|% "lane1-A2-B_S2_L001_1.trim_bismark_bt2_pe.deduplicated.wig"
out_file = argv[2]  %|%  str_replace(in_file, ".wig", "rds")
n_row    = argv[3] %>% as.integer() %|% 3400L # Number of rows of data in the output figure

### Prepare the input ####

raw_input = read_delim(in_file,
  skip = 2,  col_names = c('position', 'density')) 
missed_reads=problems() # There will be errors in raw_input (that's okay); this helps handle them
# Okay, so every N thousand rows or so, there's the line "variableStep"
# This indicates a new chromosome is about to happen

first_chrom = read_lines(in_file,  skip = 1, n_max  = 1)
# Get the chromosome details (from the parsing errors), 
# making space for the first chromosome (which was skipped)
# We're combining it with the missed reads and 
# assigning it to "row 0", since the first row in the data is actually row 2
chrom_info = tibble(row = 0, actual = first_chrom) %>% 
  bind_rows(missed_reads %>% 
            filter(expected == "a double") %>% 
            select(row, actual)) %>% 
  # Split the data into chromosome number and span
  separate(actual, into =c(NA, "chrom", "span"), 
           sep = " ") %>% 
  mutate(chrom = str_remove(chrom, fixed("chrom=")),
         span = str_remove(span, fixed("span=")) %>% as.integer())

# Add a null row to the beginning of the input data for "row 0"
count_data = tibble(position = NA_real_, density = NA_real_, row = 0) %>% 
  bind_rows(raw_input %>% mutate(row = 1:n())) %>% 
  # Then add row numbers and join against chrom_info
  # The join line will be NA for the position and density
  left_join(chrom_info, by = 'row') %>% 
  # Fill it down, then filter out the NA's
  fill(chrom, span) %>% 
  filter(!is.na(density)) %>% 
  select(-row) %>% 
# There are a number of rows with span < 25; these are the tail ends of chromosomes
# Filtering them out will be a minuscule decrease in data
  filter(span == 25) 

# Wrap the data into rows for easier plotting
out_data = count_data %>% 
  ungroup() %>%
  # select(-column) %>%
  mutate(column = rep_along(1:n(), 1:n_row) %>% sort()) %>% 
  group_by(column) %>% 
  mutate(y = 1:n())


write_rds(out_data, out_file)




