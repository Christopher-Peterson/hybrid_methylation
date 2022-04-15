
# WGBS Alignment Experiments

This script assumes that you have You have cloned the git repository,
[configured Singularity](../docker/), and set up the [WGBS data &
genome](../wgbs_setup/). Before beginning, make sure everything in
`scripts` is executable:

``` bash
cds hybrid_methylation
chmod +x alignment_experiments/scripts/*
```

## Determining Alignment Parameters

Aligning the WGBS sequencing data will be challenging, and there are a
variety of ways that it can be done. This will begin by testing various
bismark parameters.

Key metrics:

-   Mapping efficiency
-   Methylation context percentages
-   Conversion efficiency based on lambda DNA spike

First, we’ll define this function to set up a new experiment

``` bash
source alignment_experiments/functions.sh
# new_alignment_experiment name "pe" "bismark_parms"

# Alignment experiments w/ default parameters
new_alignment_experiment pe_default "pe" " "
new_alignment_experiment se_default "se" " "
new_alignment_experiment pe_relaxed "pe" "-N 1 --score_min L,0,-0.6"
new_alignment_experiment se_relaxed "se" "-N 1 --score_min L,0,-0.6"
new_alignment_experim ent pe_xrelaxed "pe" "-N 1 --score_min L,0,-0.8"
new_alignment_experiment se_xrelaxed "se" "-N 1 --score_min L,0,-0.8"
new_alignment_experiment pe_local "pe" "-N 1 --local"
new_alignment_experiment se_local "se" "-N 1 --local"
new_alignment_experiment pe_local_relax "pe" "-N 1 --local --score_min G,0,-0.6"
new_alignment_experiment se_local_relax "se" "-N 1 --local --score_min G,0,-0.6"

# Run the alignment trials
run_alignment trial pe_default
run_alignment trial se_default
run_alignment trial pe_relaxed
run_alignment trial se_relaxed
run_alignment trial pe_xrelaxed
run_alignment trial se_xrelaxed
run_alignment trial pe_local
run_alignment trial se_local
run_alignment trial pe_local_relax
run_alignment trial se_local_relax
```

### Experimental Parameters and results

The trial set for each of these is the largest trim file,
`lane1-A2-B_S2_L001_1.trim` (and its partner, for PE).

| Directory      | Single / Paired Ends? | Parameters                      | Time    | Mapping Efficiency | No Alignment | CpG Context | CHG Context | CHH Context | Unknown Context | Conversion Efficiency (lambda) | Comments                                          |
|:---------------|:----------------------|:--------------------------------|:--------|:-------------------|:-------------|:------------|:------------|:------------|:----------------|:-------------------------------|:--------------------------------------------------|
| se_relaxed     | se                    | -N 1 –score_min L,0,-0.6        | 3:22:22 | 42.4%              | 44.8%        | 9.2%        | 1.7%        | 1.9%        | 21.8%           |                                |                                                   |
| pe_relaxed     | pe                    | -N 1 –score_min L,0,-0.6        | 3:44:53 | 36.2%              | 53.7%        | 9.5%        | 1.7%        | 1.9%        | 23.6%           |                                |                                                   |
| se_default     | se                    | \[Default parameters\]          | 1:2:1   | 19.6%              | 70.0%        | 9.1%        | 0.9%        | 1.0%        | 18.5%           |                                |                                                   |
| pe_default     | pe                    | \[Default parameters\]          | 1:24:7  | 14.8%              | 77.7%        | 9.1%        | 1.0%        | 1.0%        | 19.8%           |                                |                                                   |
| se_xrelaxed    | se                    | -N 1 –score_min L,0,-0.8        | 2:42:47 | 50.3%              | 35.4%        | 9.2%        | 1.9%        | 2.2%        | 19.7%           |                                |                                                   |
| pe_xrelaxed    | pe                    | -N 1 –score_min L,0,-0.8        | 4:11:20 | 44.4%              | 44.35%       | 9.4%        | 2.0%        | 2.2%        | 21.6%           |                                |                                                   |
| se_local       | se                    | -N 1 –local                     | timeout |                    |              |             |             |             |                 |                                |                                                   |
| pe_local       | pe                    | -N 1 –local                     | 5:26:54 | 63.7%              | 25.7%        | 8.2%        | 1.0%        | 1.5%        | 4.9%            |                                |                                                   |
| se_local_relax | se                    | -N 1 –local –score_min G,0,-0.6 | 0:58:25 | 0%                 | 100%         | NA          | NA          | NA          | NA              |                                | I should figure out how score_min works for local |
| pe_local_relax | pe                    | -N 1 –local –score_min G,0,-0.6 | 0:42:36 | 0%                 | 100%         | NA          | NA          | NA          | NA              |                                |                                                   |

Based on these results, `pe_xrelaxed` is the best option that doesn’t
require `local` (e.g., soft-clipping). I’m hesitant about using that
with WGBS data.

Run the alignment on all of samples with `pe_xrelaxed`:

``` bash
# First, modify the run slurm file to have 5 hours, just to be safe
# Then, run the full alignment
run_alignment run pe_xrelaxed
```

## De-duplicate the BAM files

De-duplicate the bam files. We’ll be moving to the `pe_xrelaxed`
directory and staying there until otherwise noted.

``` bash
cds hybrid_methylation/alignment_experiments/pe_xrelaxed

sbatch slurm/bismark_dedup_run.slurm
```

Let’s verify that file sizes have changed:

``` bash
du -h bams/*bam # Original bams
# Results:
15G     bams/lane1-A1-A_S1_L001_1.trim_bismark_bt2_pe.bam
22G     bams/lane1-A2-B_S2_L001_1.trim_bismark_bt2_pe.bam
13G     bams/lane1-A3-C_S3_L001_1.trim_bismark_bt2_pe.bam
17G     bams/lane1-A4-D_S4_L001_1.trim_bismark_bt2_pe.bam
16G     bams/lane1-X1-E_S5_L001_1.trim_bismark_bt2_pe.bam
16G     bams/lane1-X3-F_S6_L001_1.trim_bismark_bt2_pe.bam
12G     bams/lane2-X2-E_S7_L002_1.trim_bismark_bt2_pe.bam
14G     bams/lane2-X4-F_S8_L002_1.trim_bismark_bt2_pe.bam
20G     bams/lane2-X5-A_S9_L002_1.trim_bismark_bt2_pe.bam
16G     bams/lane2-X6-B_S10_L002_1.trim_bismark_bt2_pe.bam
12G     bams/lane2-X7-C_S11_L002_1.trim_bismark_bt2_pe.bam
18G     bams/lane2-X8-D_S12_L002_1.trim_bismark_bt2_pe.bam

du -h bams/dedup/*bam  # de-duplicated
# Results:
11G     bams/dedup/lane1-A1-A_S1_L001_1.trim_bismark_bt2_pe.deduplicated.bam
16G     bams/dedup/lane1-A2-B_S2_L001_1.trim_bismark_bt2_pe.deduplicated.bam
9.1G    bams/dedup/lane1-A3-C_S3_L001_1.trim_bismark_bt2_pe.deduplicated.bam
13G     bams/dedup/lane1-A4-D_S4_L001_1.trim_bismark_bt2_pe.deduplicated.bam
12G     bams/dedup/lane1-X1-E_S5_L001_1.trim_bismark_bt2_pe.deduplicated.bam
11G     bams/dedup/lane1-X3-F_S6_L001_1.trim_bismark_bt2_pe.deduplicated.bam
7.8G    bams/dedup/lane2-X2-E_S7_L002_1.trim_bismark_bt2_pe.deduplicated.bam
9.6G    bams/dedup/lane2-X4-F_S8_L002_1.trim_bismark_bt2_pe.deduplicated.bam
13G     bams/dedup/lane2-X5-A_S9_L002_1.trim_bismark_bt2_pe.deduplicated.bam
11G     bams/dedup/lane2-X6-B_S10_L002_1.trim_bismark_bt2_pe.deduplicated.bam
8.5G    bams/dedup/lane2-X7-C_S11_L002_1.trim_bismark_bt2_pe.deduplicated.bam
11G     bams/dedup/lane2-X8-D_S12_L002_1.trim_bismark_bt2_pe.deduplicated.bam
```

### Extract Methylation Calls

I tested the methylation extraction script on a single de-duplicated BAM
file (`lane1-A2-B_S2_L001_1.trim_bismark_bt2_pe.deduplicated.bam`), then
verify that the output is present.

``` bash
# Run extraction trial
sbatch slurm/bismark_extract_trial.slurm
# Wait until complete
du -h methyl_extract/*lane1-A2-B_S2_L001_1.trim_bismark_bt2_pe.deduplicated*
# Results:
25G     methyl_extract/CpG_context_lane1-A2-B_S2_L001_1.trim_bismark_bt2_pe.deduplicated.txt
152M    methyl_extract/lane1-A2-B_S2_L001_1.trim_bismark_bt2_pe.deduplicated.bedGraph.gz
159M    methyl_extract/lane1-A2-B_S2_L001_1.trim_bismark_bt2_pe.deduplicated.bismark.cov.gz
957M    methyl_extract/lane1-A2-B_S2_L001_1.trim_bismark_bt2_pe.deduplicated.CpG_report.txt
2.0K    methyl_extract/lane1-A2-B_S2_L001_1.trim_bismark_bt2_pe.deduplicated.cytosine_context_summary.txt
27K     methyl_extract/lane1-A2-B_S2_L001_1.trim_bismark_bt2_pe.deduplicated.M-bias.txt
1.0K    methyl_extract/lane1-A2-B_S2_L001_1.trim_bismark_bt2_pe.deduplicated_splitting_report.txt
109G    methyl_extract/Non_CpG_context_lane1-A2-B_S2_L001_1.trim_bismark_bt2_pe.deduplicated.txt

# This is more or less what's expected, so run the rest

sbatch slurm/bismark_extract_run.slurm
```

### Compress & Backup to `WORK`

This has generated a lot of files, so it’s worth backing up the
important ones to `WORK` in case `SCRATCH` is purged.

``` bash
WRK=${PWD/$SCRATCH/$WORK} # Work version of pe_xrelaxed
mkdir -p $WRK
mkdir -p $WRK/bams/dedup

# Some files need compressing
sbatch slurm/compress_bismark_results.slurm

# Backup reports
cp -r reports $WRK/reports

# Wait on these next two until the slurm job is complete:

# Backup deduplicated BAMS
cp -r bams/dedup/*bam $WRK/bams/dedup/ &

# Backup extracted methylation information
cp -r methyl_extract/*{.gz,_summary.txt,M-bias.txt,_splitting_report.txt} $WRK/methyl_extract
```

## Sanity Checking the alignment

Before proceeding, we need to check the alignments to make sure that
they make sense. We’ll be using [IGV](https://igv.org) for this. Before
using IGV, we need to sort and index the deduplicated bam files with
samtools. (Also, the reference genome should be intexed; Go back and add
this to the genome prep script!)

### Sorting the BAM files

First, we’re using the sort & index script on a single bam to ensure the
script works.

``` bash
sbatch slurm/sort_trial.slurm
```

That ran in about 10 minutes, so let’s submit the full run.

``` bash
sbatch slurm/sort_run.slurm
```

### Installing IGV

First, download and unzip IGV to TACC.

``` bash
cdw 
wget https://data.broadinstitute.org/igv/projects/downloads/2.12/IGV_Linux_2.12.3_WithJava.zip 
unzip IGV_Linux_2.12.3_WithJava.zip
rm IGV*zip
```

Next, we want to make sure that IGV doesn’t start writing data to $HOME
when we run it. We’ll do this by making a scratch directory for IGV and
then linking it to home.

``` bash
cdw
mkdir IGV_genomes
cdh
mkdir igv
cd igv
ln -s $WORK/IGV_files genomes
```

### Using IGV

We’re going to use the [TACC Visualization Portal](vis.tacc.utexas.edu/)
to use IGV. Go to the portal and start a DCV remote desktop session on
Lonestar 6 using the development queue. If DCV is unavailable, you can
use VNC instead. When the page changes, click the green “Connect” button
(once it appears), then sign into TACC on the DCV page.

On the virtual desktop, open a terminal, navigate to where you saved
IGV, and run `igv.sh`. If you want, you can right-click on the desktop
and create a launcher (a.k.a., a shortcut) for it to save time in the
future.

# ToDo:

Sanity checks with IGV viewer Lambda Spike?

SNPsplit
