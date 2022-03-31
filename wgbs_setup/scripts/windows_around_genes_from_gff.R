#!/usr/bin/env Rscript
#windows_around_genes_from_gff.R
#This builds a bed file indicating windows upstream and downstream from genes
#output is a bed file with first three columns as required,
#and a final 'name' column with information about each window.
#In the names, windows are numbered relative to the GENE'S start

#Example names for a forward gene:
#w1-Amillepora16770-574-1574-upstream
#w79-Amillepora16770-78574-79574-inside
#w86-Amillepora16770-86389-87389-downstream

#Example names for a forward gene:
# w1-Amillepora16774-268672-269672-upstream
# w101-Amillepora16774-168672-169672-inside
# w141-Amillepora16774-128293-129293-downstream

#PARSE ARUGMENTS
suppressPackageStartupMessages({
  library(tidyverse)
  library(optparse)
  library(parallel) 
  })

option_list = list(
  
  make_option(c("--i"), type="character", default=NULL, 
              help="infile"),
  
  make_option(c("--gene_id"), type="character", default='ID', 
              help="gene_id"),
  
  make_option(c("--window_size"), type="integer", default=1000, 
              help="size of the windows in bp"),
  
  make_option(c("--buffer"), type="integer", default=50000, 
              help="how far out upstream and downstream to make windows from gene starts"),
  
  make_option(c("--o"), type="character", default=NULL, 
              help="output name"),
  make_option(c("--cores"), type = "integer", default = 1, 
              help = "number of cores to use")
)

print("Parsing arugments...")
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

infile = opt$i
gene_id = opt$gene_id
window_size = opt$window_size
buffer = opt$buffer
outname = opt$o
N_CORES = opt$cores

print('Reading in GFF...')
gff_raw = read_tsv(infile,comment = "#",
               col_names = c('chr', 'source', 'feature', 'start', 'end', '.', 'strand', '..', 'description'))

print('Getting maximum lengths for chromosomes...')
#get maximum lengths
max_length_df = gff_raw %>% 
  group_by(chr) %>% 
  summarize(max_length = max(end))
max_length = max_length_df$max_length
names(max_length) = max_length_df$chr

#filter to only genes
gff = gff_raw %>% 
  filter(feature=='gene') %>% 
  mutate( # Convert these types so they work with the functions below
    start = as.numeric(start),
    end = as.numeric(end),
    chr = as.character(chr),
    description = as.character(description),
    gene = description %>% str_remove(paste0("^.*", gene_id, "\\=")) %>% str_remove("\\;.*$"),
    max_len = max_length[chr],
    buffer_left = pmax(start - buffer, 0),
    buffer_right = pmax(end + buffer, 0)
  ) %>% select(-description, -`.`, -`..`)

# HANDLE THE FORWARD GENES ------------------------------------------------

#if the set of boundaries is too short, just return NULL
window_from_bounds = function(bounds, region, chr, gene, direction = c("forward", "reverse")){
  # If direction is rev, switch r and l
  if (length(bounds) == 0) return(NULL)
  fwd_l = bounds[1:(length(bounds)-1)]
  fwd_r = bounds[2:length(bounds)]
  name = paste(gene, region, sep='_')
  if(direction == 'forward') {
    w = tibble(chr, l=fwd_l, r = fwd_r, name = name)
  } else {
    w = tibble(chr, l=fwd_r, r = fwd_l, name = name)
  }
  w
}

# Faster version that can be parallelized
forward_gene_window_row = function(row){
  # row should have columns start (numeric), end(numeric), and chr
  with(row, {
    left_bounds = rev(seq(start, buffer_left, by = -window_size))
    inside_bounds = seq(start, end, by = window_size)
    inside_bounds[length(inside_bounds)] = end #add leftover into final window
    right_bounds = seq(end, buffer_right, by = window_size)
    
    left_windows = window_from_bounds(left_bounds, 'upstream',chr, gene, direction = 'forward')
    inside_windows = window_from_bounds(inside_bounds, 'inside',chr, gene, direction = 'forward')
    right_windows = window_from_bounds(right_bounds, 'downstream',chr, gene, direction = 'forward')
    
    # dist = res0$l - start # I don't think this does anything
    tss = paste('tss',start, sep='')
    tts = paste('tts', end, sep='')
    
    bind_rows(left_windows, inside_windows, right_windows) %>% 
      mutate(
        window_nums = paste('w', 1:n(), sep=''),
        name = paste(name, tss, tts, sep = "_")
    ) # return
  })
}
reverse_gene_window_row = function(row){
  # row should have columns start (numeric), end(numeric), and chr
  with(row, {
    left_bounds = seq(start, buffer_left, by = -window_size) #NOTE DIFFERENCE FROM forward_gene_windows
    inside_bounds = seq(end, start, by = -window_size) #NOTE DIFFERENCE FROM forward_gene_windows
    inside_bounds[length(inside_bounds)] = start #include any leftover
    right_bounds = rev(seq(end, buffer_right, by = window_size)) #NOTE DIFFERENCE FROM forward_gene_windows
    
    left_windows = window_from_bounds(left_bounds, 'downstream',chr, gene, direction = 'reverse')
    inside_windows = window_from_bounds(inside_bounds, 'inside',chr, gene, direction = 'reverse')
    right_windows = window_from_bounds(right_bounds, 'upstream',chr, gene, direction = 'reverse')

    tss = paste('tss',start, sep='')
    tts = paste('tts', end, sep='')
    
    bind_rows(left_windows, inside_windows, right_windows) %>% 
      mutate(
        window_nums = paste('w', 1:n(), sep=''),
        name = paste(name, tss, tts, sep = "_")
      ) # return
  })
}


print('Buiding windows for forward genes...')
forward_bed = gff %>% 
  filter(strand == '+') %>% 
  rowwise() %>% group_split() %>% # Split into a list of 1-row df's
  mclapply(forward_gene_window_row, mc.cores = N_CORES) %>% 
  bind_rows()

print('Buiding windows for reverse genes...')
reverse_bed = gff %>% 
  filter(strand == '-') %>% 
  rowwise() %>% group_split() %>% # Split into a list of 1-row df's
  mclapply(reverse_gene_window_row, mc.cores = N_CORES) %>% 
  bind_rows()

# output ------------------------------------------------------------------

print(paste('Saving results as ', outname, '...', sep=''))
bed = bind_rows(forward_bed, reverse_bed)
bed %>% 
  write_tsv(outname, col_names = FALSE)
