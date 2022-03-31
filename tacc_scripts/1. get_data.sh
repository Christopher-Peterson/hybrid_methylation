# Hybrid Methylation setup

# In this doc:
  # Install the singularity images
  # Setup the reference genome
  # Get the data, unpack it, organize it
  # create a small subset of the fasta files for QC testing
  # Trim the fasta files
  # QC test the trimmed ones

# First, make sure that the basic files have been copied over

sbatch slurm/setup_singularity.sh
# Wait until this one is done, then run the reference genome setup
sbatch slurm/setup_reference_genome.sh

source singularity_aliases.sh

#### Acquire the files ######
# Copy data from Grove's archive
cds hybrid_methylation/raw_data
cp  /corral-repl/utexas/tagmap/dixon_backups/file_sharing_test_dir/JA20120.tar.gz .

# Untar
tar –xvzf JA20120.tar.gz
cd JA20120_download
# There should be 12 directories in here
ls | wc -l # 13
rm *json # Remove an unneeded json file

# Each of these directories containes two fasta.gz files (the paired ends)
# Move them to a new, single directory
mkdir ../fastqs
mv */*gz ../fastqs
cd ../fastqs
ls
ls | wc -l # should be 24

# Now, I should unzip all of these
gunzip *.gz
# This was taking too long, so I opened up idev and tried it in parallel;
  # the first three files were already done
idev
pigz -d lane1-A2-B_S2_L001_R2_001.fastq.gz &
pigz -d lane1-A3-C_S3_L001_R1_001.fastq.gz &
pigz -d lane1-A3-C_S3_L001_R2_001.fastq.gz &
pigz -d lane1-A4-D_S4_L001_R1_001.fastq.gz &
pigz -d lane1-A4-D_S4_L001_R2_001.fastq.gz &
pigz -d lane1-X1-E_S5_L001_R1_001.fastq.gz &
pigz -d lane1-X1-E_S5_L001_R2_001.fastq.gz &
pigz -d lane1-X3-F_S6_L001_R1_001.fastq.gz &
pigz -d lane1-X3-F_S6_L001_R2_001.fastq.gz &
pigz -d lane2-X2-E_S7_L002_R1_001.fastq.gz &
pigz -d lane2-X2-E_S7_L002_R2_001.fastq.gz &
pigz -d lane2-X4-F_S8_L002_R1_001.fastq.gz &
pigz -d lane2-X4-F_S8_L002_R2_001.fastq.gz &
pigz -d lane2-X5-A_S9_L002_R1_001.fastq.gz &
pigz -d lane2-X5-A_S9_L002_R2_001.fastq.gz &
pigz -d lane2-X6-B_S10_L002_R1_001.fastq.gz &
pigz -d lane2-X6-B_S10_L002_R2_001.fastq.gz &
pigz -d lane2-X7-C_S11_L002_R1_001.fastq.gz &
pigz -d lane2-X7-C_S11_L002_R2_001.fastq.gz &
pigz -d lane2-X8-D_S12_L002_R1_001.fastq.gz &
pigz -d lane2-X8-D_S12_L002_R2_001.fastq.gz &
# monitor with top
# exit idev w/ exit

# Check numbers
ls *.fastq | wc -l	# 24
ls lane1*.fastq | wc -l # 12
ls lane2*.fastq | wc -l # 12
ls *_R1_*.fastq | wc -l # 12
ls *_R2_*.fastq | wc -l # 12



cds hybrid_methylation
mkdir reads
cd reads
ln -s ../raw_data/fastqs raw
cd ..

