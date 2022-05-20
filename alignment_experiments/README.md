
# WGBS Alignment Experiments

This script assumes that you have You have cloned the git repository,
[configured Singularity](../setup/docker/), set up [IGV](../setup/igv),
and prepared the [WGBS data & genome](../setup/wgbs_setup/). Before
beginning, make sure everything in `scripts` is executable:

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
with WGBS data. However, additional diagnostics are required before
running everything.

### De-duplicate the BAM files

We need to de-duplicate the bam files from the trials that seemed
vaguely promising.

``` bash
dedup_alignment pe_xrelaxed
dedup_alignment pe_relaxed
dedup_alignment se_relaxed
dedup_alignment se_xrelaxed
# dedup_alignment pe_local
```

Let’s verify that file sizes have changed in one of the directories:

``` bash
cds hybrid_methylation/alignment_experiments/pe_xrelaxed

du -h bams/*bam # Original bams
# Results:
22G     bams/lane1-A2-B_S2_L001_1.trim_bismark_bt2_pe.bam

du -h bams/dedup/*bam  # de-duplicated
# Results:
16G     bams/dedup/lane1-A2-B_S2_L001_1.trim_bismark_bt2_pe.deduplicated.bam
cd ..
```

### Create the tapestry plots

Tapestry plots are a visualization method for checking how read coverage
varies across the genome. This process involves:

-   Sorting and indexing the deduplicated BAM files using samtools.  
-   Estimating read counts across the genome with IGVtools.
-   Visualizing the read counts with R and ggplot2.

(Also, the reference genome should be intexed; Go back and add this to
the genome prep script!)

``` bash
vis_tapestry pe_xrelaxed
vis_tapestry pe_relaxed
vis_tapestry se_relaxed
vis_tapestry se_xrelaxed
# vis_tapestry pe_local
```

# ToDo:

-   Lambda Spike?
-   SNPsplit

## Run full alignment (have been using pe_xrelaxed, but may need to switch)

Run the alignment on all of samples with `pe_xrelaxed`:

``` bash
# First, modify the run slurm file to have 5 hours, just to be safe
# Then, run the full alignment
run_alignment run pe_xrelaxed
```

We’ll be moving to the `pe_xrelaxed` directory and staying there until
otherwise noted. ### Extract Methylation Calls

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
