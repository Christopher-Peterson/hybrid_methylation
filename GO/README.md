
# Gene Ongology analysis

The goal here is to see GO terms shed light on the mixed heritability
(as revealed in DSS).

## Get Gene Names of DML

``` r
# From the full dss output bed file
dat |> select(1:4) |> distinct() |>  write_tsv('GO/dml_pos.bed')
```

Push it up to TACC, then attach gene names

``` bash

window_dir=$STOCKYARD/tagmap-share/genomes/Amil_v2.1_ForBismark
db="${window_dir}/geneBoundaries.bed"

echo -e '#chrom\tstart\tstop\toffspring\t#chrom2\tstart2\tstop2\tgene' > dml_genes.bed
bedtools intersect -a dml_pos.bed -b $db -wa -wb >> dml_genes.bed
head dml_genes.bed
```

Pull it into R to examine it

``` r
dml_genes = read_tsv('GO/dml_genes.bed')
# dml_genes |> filter(`#chrom` != `#chrom2`) # sanity check
dml_trim = dml_genes |> select(1:4, gene_id = gene) # Other starts are for genome feature

# Average heritability per-gene
dss_data = read_tsv('dss/model_data/dss_filtered_uncorrelated_0.25_N8.bed')
gene_joint = dml_trim |> 
  left_join(dss_data |> select(-parent, -role, -delta_po, -delta.se_po) |> distinct())
gene_herit = gene_joint |> mutate(eta = diff_delta * sign(delta_p)) |> 
  group_by(gene_id) |> 
  summarize(mean_eta = mean(eta), mean_abs_eta = mean(abs(eta)),
            sd_eta = sd(eta), sd_abs_eta = sd(abs(eta)), n = n()) |> 
  ungroup() |> 
  mutate(gene_id = str_remove(gene_id, '^gene-'))

write_csv(gene_herit, 'GO/db/gene_herit.csv')
# This is the version that will actually be used.
gene_herit |> 
  select(gene = gene_id, eta = mean_eta) |> 
  mutate(eta = abs(eta)) |> 
  write_csv('GO/db/gene_herit_mean.csv')
# gene_herit |> ggplot(aes(x = gene_id, y = mean_eta)) + 
#   geom_point()

# Alternate version of heritability:
  # number of loci in the gene that were ID'd as heritable vs. non-heritable in the locus-specific analysis
```

Alternately, see if we can skip steps 1-2 by downloading the ensembl
file for Amil v 2.1.

Drop all lines that don’t contain: ‘join(’, ‘/gene’, ‘/db_xref=“GO:’

``` bash
mkdir -p GO/db
cd GO/db
wget https://ftp.ensembl.org/pub/rapid-release/species/Acropora_millepora/GCA_013753865.1/refseq/geneset/2022_03/Acropora_millepora-GCA_013753865.1-2022_03-genes.embl.gz

idev
gunzip *gz


extract_go_genes() {
  # Create a TSV w/ gene names & associated GO terms
  # Filter for gene names, db refs, and GO terms
  local filter_rows='join\(|/gene|/db_xref.\"GO'
  # Trim headers and quotes
  local trim1='s%FT */%/%' # must be in quotes
  local trim2='s/\"//g'
  # combine each record onto one line
  #local merge_lines1="tr --delete '\n'" # commented out because it didn't work in quotes for some reason
  local merge_lines2='s/FT/\n/g'
  # Subset to only included the gene immediately before GO term(s) (which correspond to them)
  local get_go='/gene\=LOC[0-9]+(/db_xref\=GO\:[0-9]+)+' 
  # format as tsv, in four steps
  # Remove the first gene=
  local format_tsv1='s%/gene=%%'
  # Newline for the remaining gene=
  local format_tsv2='s%/gene=%\n%g'
  # Tab before the first go term
  local format_tsv3='s%/db_xref=%\t%'
  # semicolons subsequently
  local format_tsv4='s%/db_xref=%;%g'
  
 grep -P "$filter_rows" | \
    sed -e "${trim1}" -e "${trim2}"  | tr --delete '\n' | \
    sed -e "${merge_lines2}" | \
    grep -oP $get_go | sed -e "$format_tsv1" -e "$format_tsv2" | \
    sed -e "$format_tsv3" -e "$format_tsv4"
}; 

# These are the headers compatible w/ GO-MWU
echo -e 'V1\tV2' > Amil_2.1_GO.tsv
cat Acropora_millepora-GCA_013753865.1-2022_03-genes.embl | extract_go_genes >> Amil_2.1_GO.tsv
cd ..

# Columns are Gene_ID \t comma-separated_list_of_GO_terms
```

I’ve downloaded and unzipped the GO_MWU repository into the GO
directory; some adjustments have been made to the script re: filenames.

Acquire the go.obo database

``` bash
cd GO/db
wget https://purl.obolibrary.org/obo/go.obo
cd ../GO_MWU-master

# Link the files from db
ln -s ../db/Amil_2.1_GO.tsv Amil_GO.tsv
ln -s ../db/gene_herit_abs.csv gene_herit_abs.csv
ln -s ../db/gene_herit_mean.csv gene_herit_mean.csv

ln -s ../db/go.obo go.obo

Rscript GO_MWU.R gene_herit_abs.csv Amil_GO.tsv BP
```

``` r
# Use withr::with_wd() to run gomwu?  Or do it cmd.
# yeah, that
```

<!-- Next step: use bedtools  -->
<!-- Run dss/scripts/dml_rds_to_bed.r, which converts each chrom * pair rds file into a bed file.  -->
