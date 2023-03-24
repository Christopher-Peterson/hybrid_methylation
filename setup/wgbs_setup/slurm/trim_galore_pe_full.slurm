#!/bin/bash

#-------------------------------------------------------
#SBATCH -J trim_galore_pe
#SBATCH -N 4 
#SBATCH -n 12
#SBATCH -p development
#SBATCH -o logs/trim_galore_pe.o
#SBATCH -e logs/trim_galore_pe.e
#SBATCH -t 00:40:00
### #SBATCH -A tagmap # Seems like we're not using tagmap any more

#------------------------------------------------------
# TRIM PAIRED-END READS

module load tacc-singularity
module unload xalt
module load launcher
export LAUNCHER_WORKDIR=`pwd`

CORES_PER_JOB=42
# Populate the job file
export LAUNCHER_JOB_FILE=jobs/trim_galore_pe_full
: > $LAUNCHER_JOB_FILE

in_dir=raw_data/fastqs
out_dir=raw_data/trimmed
qc_dir=raw_data/fastqc
mkdir -p  $out_dir $qc_dir

function add_job {
  file_2=$1
  file_1=${file_2/_R2_/_R1_}
  # out_1=${file_2/_R2_001.fastq/}_1.trim
  # out_2=${file_2/_R2_001.fastq/}_2.trim
  # The clip_R1 options are recommended for picomethyl
  echo "trim_galore -q 20 --fastqc \
  --fastqc_args \"--outdir ${qc_dir}\" \
  --clip_R1 10 --clip_R2 10 --three_prime_clip_R1 10 \
  --three_prime_clip_R2 10 \
  -j $CORES_PER_JOB --paired \
  -o $out_dir $file_1 $file_2 " >> $LAUNCHER_JOB_FILE
}

# Create the job script
for fq_r2 in $in_dir/*_R2_001.fastq; do
  add_job $fq_r2
done

# Run the launcher
$LAUNCHER_DIR/paramrun

# echo "Done with run"
# cds hybrid_methylation/reads
# mkdir trimmed
# mv raw/*trim trimmed
# rm trimmed/*trim.temp.*
