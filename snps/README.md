
# SNP extraction from alignment

This assumes the alignment pe_xrelaxed has been completed in the
[Alignment Experiments section](../alignment_experiments).

Essentially, this task aims to:

1.  Mask the alignment so it can have variants called (w/ `revelio.py` )
2.  Call SNPs from the alignment files (making VCF files)
3.  Create N-masked genomes from the VCF files, for use in SNPsplit

## Setup data & files

First, let’s set up the file structure.

``` bash
cd $SCRATCH/hybrid_methylation
mkdir snps snps/scripts snps/slurm
cd snps/scripts
# Download the masking script
wget https://raw.githubusercontent.com/bio15anu/revelio/main/revelio.py

cd ..
```

Setup the folder for the alignment.

``` bash
ALIGNMENT=pe_xrelaxed
mkdir $ALIGNMENT
cd $ALIGNMENT
mkdir jobs split_genomes temp vcf logs masked_bams masked_genomes offspring offspring/bams offspring/reports offspring/snpsplit offspring/snps
chmod +x scripts/SNPsplit_genome_preparation
chmod +x scripts/*.sh
# Link the sorted bams
ln -sf $SCRATCH/hybrid_methylation/alignment_experiments/$ALIGNMENT/bams/sorted bams
# Copy the genome link
ln -sfL $STOCKYARD/hybrid_methylation/alignment_experiments/$ALIGNMENT/genome.sh genome.sh
# Link the config from alignment
ln -sf $SCRATCH/hybrid_methylation/alignment_experiments/$ALIGNMENT/config config


# Link the slurm & script folders
ln -s ../slurm slurm
ln -s ../scripts scripts
```

## Generate VCF files for parents

Now, we’ll generate 6 copies of the genome (one for each parent
combination) where fixed SNP differences between the parents are coded
as N (ambiguous); the offspring reads will be re-mapped to the masked
genomes of their parents and SNP-split can use this to figure out
allele-specific methylation.

First, we need to generate VCF files for each pair of parents. We use a
pre-processing step that reduces C-\>T conversion bias inherent in in
BS-Seq data (revelio.py); then, we use ANGSD to create and merge the VCF
files. Finally, we filter the VCFs to only include SNPS where both
parents had fixed, opposite genotypes and genotype probabilities were \>
0.95.

``` bash
sbatch slurm/revelio_run.slurm # Creates masked bam files`

# These masked Bams need to be indexed
# FOr whatever reason, the revelio script didn't do the indexing
for bam in masked_bams/*bam; do
samtools index $bam &
done

# Once this is done, run this sequence of slurm scripts
function sbatch_job {
  local job=`sbatch "${@}" | grep "Submitted batch job" | grep -Po "[0-9]+$" `
  echo $job
}

# Run the VCF script
j0=`sbatch_job -d  slurm/angsd_vcf_run.slurm`
# Merge the results
j1=`sbatch_job -d afterok:$j0 slurm/angsd_vcf_merge.slurm`
# Filter the VCFs to only include snps w/ fixed, opposite genotypes & two samples
j2=`sbatch_job -d afterok:$j1 slurm/vcf_filter.slurm`
```

# Combining mdRAD and WGBS data to improve SNP discrimination

``` bash
cds hybrid_methylation/snps

mkdir -p with_md with_md/jobs with_md/logs
cd with_md
ALIGNMENT=pe_xrelaxed
# Link appropriate files & directories
ln -sf ../$ALIGNMENT/masked_bams wgbs_masked_bams
ln -sf ../../setup/md*/bams md_bams
ln -sf ../slurm slurm
ln -sf ../scripts scripts
ln -Lsf ../$ALIGNMENT/genome.sh genome.sh
ln -sf ../../$ALIGNMENT/vcf/concat vcf/wgbs
```

Note, the naming conventions between the WGBS and mdRAD data are
different; here’s the adult conversion:

- A1 = selago4
- A2 = selago6
- A3 = milepora8
- A4 = milepora11

And here’s the association between parents and offspring:

``` bash
# Offspring, Parent Combo
echo "Offspring Parent
X1    A12   
X3    A34   
X5    A13   
X6    A14   
X7    A23   
X8    A24   
X2    A12   
X4    A34" > parents.txt
```

Now, we’re going to call SNPs for the mdRAD data and merge the SNP data
into the WGBS VCF.

``` bash
# Offspring, Parent Combo

mkdir -p offspring offspring/bams offspring/reports offspring/snpsplit offspring/snps

# Set up the config params, but with some modifications for N-masking
ORIG_PARMS=`cat  ../$ALIGNMENT/config/bismark.parms`
mkdir -p config
echo "$ORIG_PARMS -np 0" > config/bismark.parms

# Run the VCF script on the mdRAD data

# combine the paired MD bam files
j3=`sbatch_job -d afterok:${j2} slurm/concat_parent_mds.slurm`
j4=`sbatch_job -d afterok:${j3} slurm/angsd_vcf_md.slurm`
# Merge the results
j5=`sbatch_job -d afterok:${j4} slurm/angsd_vcf_merge.slurm`

# Filter and join the WGBS & mdRAD data + mask the parental genomes
j6=`sbatch_job -d afterok:${j5} slurm/vcf_filter_join_mask.slurm`
```

Now that parental genomes are masked, we perform allele-specific
alignment with the offspring WGBS reads.

``` bash
# run new alignments 
j7=`sbatch_job -d afterok:${j6}  slurm/align_masked.slurm`
#j4=`sbatch_job slurm/align_masked.slurm`

# Dedup, sort, index alignments
j8=`sbatch_job -d afterok:${j7}  slurm/dedup_sort_offspring.slurm`

# Now run the snpsplit script
j9=`sbatch_job -d afterok:${j8}  slurm/snpsplit_trial.slurm`

# Methylation extract
j10=`sbatch_job -d afterok:${j9}  slurm/methyl_extract.slurm`
```

## Setup for differential methylation analysis

Finally, we want to identify how the order of reference/alternate genome
matches to dam/sire; this won’t be consistent between individual
offspring because we used the same reference genome for each pair of
pure crosses.

    ## stat: cannot stat 'offspring/snpsplit/X*tsv': No such file or directory

``` bash
# Link to the parents cov.gz files
mkdir -p offspring/snpsplit/methyl_extract/parents
for i in {1..4}; do
ln -sLf $SCRATCH/hybrid_methylation/alignment_experiments/pe_xrelaxed/methyl_extract/lane1-A${i}*.cov.gz offspring/snpsplit/methyl_extract/parents/A${i}.cov.gz
done
```

# Pipeline counts

``` bash
cds hyb*/snps/with_md
pc_bam=pipeline_counts_offspring_bam.tsv
pc_dedup=pipeline_counts_offspring_dedup.tsv
pc_snpsplit=pipeline_counts_offspring_snpsplit.tsv
snp_df=pipeline_snp_counts.tsv
# Main bam
bam_pat="(Sequence pairs analysed in total|Number of paired-end alignments with a unique best hit)"
echo -e "#ID\tcolumn\tcount" > $pc_bam
grep -P "$bam_pat" offspring/reports/*trim_bismark*txt  | \
  sed  -e 's|offspring/reports/lane.-|| # remove dirname' \
  -e 's|-.*.txt:|\t| # Remove tail' \
  -e 's|[[:blank:]][[:blank:]]\+|\t| # Trim excess white space' >> $pc_bam
head $pc_bam

# dedup
# Dedup report
dedup_pat='(Total number of alignments|Total count of deduplicated leftover sequences)'
echo -e "#ID\tcolumn\tcount" > $pc_dedup
grep -P "$dedup_pat" offspring/bams/dedup/*txt | \
  sed -e 's|offspring/bams/dedup/lane.-|| # remove dirname' \
  -e 's|-.\+.txt:|\t| # remove tail of file name' \
  -e 's| in .\+.bam|| # Remove bam name' \
  -e 's/:/\t/ # fix field separation' \
  -e 's/[0-9]\+.[0-9]\+% of total//' -e 's|()|| # Remove parenthetical part' \
  -e 's|[[:blank:]][[:blank:]]\+|\t| # Trim excess white space' >> $pc_dedup

# SNPSPLIT bams

split_pat='(Read pairs.singletons processed in total|specific for genome)'
echo -e "#ID\tcolumn\tcount" > $pc_snpsplit
grep -P "$split_pat" offspring/snpsplit/*sort.txt | \
  sed -e 's|offspring/snpsplit/||' -e 's|.SNPsplit_sort.txt:|\t| # trim file name' \
  -e 's/[0-9]\+.[0-9]\+%//' -e 's|()|| # Remove parenthetical part' \
  -e 's/Reads were specific for genome /genome_/' -e 's/:/\t/ # This line and the next remove extraneous text' \
  -e 's|Read pairs/singletons processed in total|total_pairs| ' \
  -e 's|[[:blank:]][[:blank:]]\+|\t| # Trim excess white space' >> $pc_snpsplit
  
# Snp counts
echo -e "#ID\tcolumn\tcount" > $snp_df
grep -P "SNPs stored in total:" offspring/snpsplit/*report.txt | \
sed -e 's|offspring/snpsplit/||' \
-e 's|.SNP.\+total:|\tSNPs\t| # Trim excess text & replace with "SNPS" column' >> $snp_df
```
