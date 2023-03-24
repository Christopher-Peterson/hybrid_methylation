#!/bin/bash

#-------------------------------------------------------
#SBATCH -J snpsplit_trial
#SBATCH -N 1 
#SBATCH -n 1
#SBATCH -p development
#SBATCH -o logs/snpsplit_trial.o
#SBATCH -e logs/snpsplit_trial.e
#SBATCH -t 02:00:00

# ml launcher 
module unload xalt

WD=$PWD
OUTDIR=$WD/offspring/snpsplit
mkdir -p $OUTDIR
cd $OUTDIR
ln -sLf $WD/scripts scripts

function run_snpsplit {
  local ID=X${1}
  local BAM=$WD/offspring/bams/lane*${ID}*.bam
  local PARENTS=`cat $WD/parents.txt | grep $ID | grep -o A..`
  # Combine SNPS
  # May ned a different option
  local SNP_TXT=$WD/masked_genomes/${PARENTS}/all_*based_on*.txt
  
  if [ ! -f ${ID}.SNPsplit_sort.txt ]; then 
    # Link input files to output directory
    ln -sf ${BAM} ${OUTDIR}/${ID}.bam
    ln -sf ${SNP_TXT} ${OUTDIR}/${ID}_SNP.tsv
    scripts/SNPsplit --bisulfite --paired --snp_file $OUTDIR/${ID}_SNP.tsv $OUTDIR/${ID}.bam
  fi  
}

for i in {1..8} ; do
  run_snpsplit $i &
done

wait
