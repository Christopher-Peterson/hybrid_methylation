########## Subset the data for QC testing #####
source singularity_aliases.sh

NTEST=1000000

mkdir ../../testrun
mkdir ../../testrun/fastqs

for file in *.fastq
do head -n $NTEST $file > ../../testrun/fastqs/${file}
done
cd ../../testrun
ls -lh fastqs

##### Run FAST QC #####
# Setup structure for slurm run
mkdir logs
mkdir jobs
mkdir slurm
mkdir fastqc_results_raw
# Upload slurm/fastqc_test.slurm to slurm/fastqc.slurm
sbatch slurm/fastqc.slurm
# Wow, that took around 2 minutes

###### Trim for Paired end reads ####

# Trim the paired-end reads
# Note: currently doing this for the test set
# I'm pretty sure it's actually supposed to happen
# on the whole datasets

# Run the slurm file
sbatch slurm/trim_pe.slurm

# This worked, in that it created output files
# Now let's try it on the real data

cd ..
mkdir trimming
cd trimming
ln -s ../raw_data/fastqs fastqs
ls
mkdir logs
mkdir jobs
mkdir slurm
# Upload over the trim_pe_full.slurm file
# This has been modified touse 42 cores per fasta pair
vim slurm/trim_pe.slurm
# :q
sbatch slurm/trim_pe.slurm

# It took less than 20 minutes to run

#### Next: Re-run QC tests on NEW subsets

# Remove old subsets and make new ones
rm ../testrun/fastqs_trimmed/*
  
  for file in fastqs_trimmed/*trim
do head -n $NTEST $file > ../testrun/${file}
done

cds hybrid_methylation/reads
ls ../trimming/fastqs_trimmed
ln -s ../trimming/fastqs_trimmed trimmed
ls trimmed
rm trimmed/*trim.temp.*
rm -r trimmed/tmp
ls
