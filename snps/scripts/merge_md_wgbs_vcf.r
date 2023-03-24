# R helper for merging md and wgbs vcf files
# Do not invoke directly; it should be called by the similarly named shell script
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})
# argv = c('/tmp/A13/md_inter.vcf', '/tmp/A13/wgbs_inter.vcf', 'vcf/concat/A13_filtered.vcf')
if(!exists('argv')) argv = commandArgs(TRUE)
md_vcf = argv[1]
wg_vcf = argv[2]
out_vcf = argv[3] # The header is already written


threads = 10
md_data = read_tsv(md_vcf,col_types = 'cccccccccc', comment = "##", num_threads = threads )
wg_data = read_tsv(wg_vcf,col_types = 'cccccccccc', comment = "##", num_threads = threads )
header = read_lines(out_vcf)

# get chrom orders
chroms = header |> str_subset("##contig") |> str_remove(".+ID\\=") |> str_remove(",.+")

# ID the lines in wg_data that aren't also in md_data
md_cols = md_data |> select(1,2) # ideal?
wg_unique = anti_join(wg_data, md_cols, by = c("#CHROM", "POS"))

# Figure out the chrom ordering from the header
#fix_sample_name = \(x) x |> str_remove('masked_bams/') |> str_remove('_masked.bam')

combined_data = bind_rows(md_data, wg_unique) |> 
  # Fix the sample names
 # rename_with(fix_sample_name, .cols = starts_with('masked_bams')) |> 
  # Convert types to allow for sorting
  mutate(`#CHROM` = `#CHROM` |> factor(levels = chroms),
         POS = as.integer(POS)) |> 
  arrange(`#CHROM`, POS)

write_tsv(combined_data, out_vcf, append = TRUE, col_names = TRUE)
