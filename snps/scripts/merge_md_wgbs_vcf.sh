#!/bin/bash

# merge and filter the bcf files from the md-RAD and WGBS data

NAME=$1 # e.g, A12
BCF_MD=$2
BCF_WG=$3

GP_MIN=.8
# NAME=A13
# BCF_MD=vcf/concat/A13.bcf
# BCF_WG=../pe_xrelaxed/vcf/concat/A13.bcf
# 
VCF_OUT=vcf/concat/${NAME}_filtered.vcf
mkdir /tmp/$NAME

N_CORE=10
## Step 1: make intermediate filtered files for each:
VCF_MD="/tmp/$NAME/md_inter.vcf"
VCF_WG="/tmp/$NAME/wgbs_inter.vcf"
# Filter the BCF files to remove heterozygotes and genotype probabilities < 0.9
bcf_filt="FMT/GP>${GP_MIN} && INFO/NS=2 && FMT/GT=\"1/1\""
# Then remove T/C and G/A SNPS 
grep_filt='T\tC|G\tA|C\tT|A\tG' # match C/T and A/G SNPS, which will be filtered out

# rename samples
P1="${NAME/%?}"
P2="A${NAME/#??}"
samp_names=/tmp/$NAME/rename
# This could be an issue...
echo "masked_bams/${P1}_masked.bam $P1
masked_bams/${P2}_masked.bam $P2
md_bam_merge/${P1}.bam $P1
md_bam_merge/${P2}.bam $P2
" > $samp_names
# cat $samp_names
# Run the intermediate filter
function filter_intermediate {
  local IN=$1
  local OUT=$2
  
  bcftools view --genotype "^het" -v snps -i "${bcf_filt}" --threads $N_CORE ${IN} | \
  grep -Pv ${grep_filt} | \
  bcftools reheader -s $samp_names > $OUT
}

filter_intermediate $BCF_MD $VCF_MD & 
filter_intermediate $BCF_WG $VCF_WG &
wait
# Write out the vcf header to VCF_OUT
grep -P "^##" $VCF_MD > $VCF_OUT
# Run the R script that will merge the two files
r-tidy-opt scripts/merge_md_wgbs_vcf.r $VCF_MD $VCF_WG $VCF_OUT

# Mask the genome
scripts/run_SNPsplit_genome_prep.sh $VCF_OUT Amil
