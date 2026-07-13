#!/bin/bash
# 01b_fixup_gzip_z01_z02.sbatch  --  finish the stalled plain-gzip on Z01/Z02 R3
# The io-array plain-gzip crawled (~0.5 MB/s); R3.fastq is fully extracted but only
# partially compressed. Remove the partial .gz and re-compress R3 with pigz (fast).
#
#SBATCH --job-name=pigzfix
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --array=1-2
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=Analysis/ATAC/Human_External/logs/fixgz_%A_%a.out
#SBATCH --error=Analysis/ATAC/Human_External/logs/fixgz_%A_%a.err
set -eo pipefail
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$PROJ"
PIGZ=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/pigz
DONORS=(Z01 Z02); D=${DONORS[$((SLURM_ARRAY_TASK_ID-1))]}
DIR="$PROJ/data/raw/GSE281367/fastq/$D"
echo "[$(date)] fixup gzip $D on $(hostname)"
cd "$DIR"
if [ -f "${D}_S1_L001_R3_001.fastq" ]; then
  rm -f "${D}_S1_L001_R3_001.fastq.gz"           # drop the partial/stalled .gz
  "$PIGZ" -p "$SLURM_CPUS_PER_TASK" "${D}_S1_L001_R3_001.fastq"
  echo "[$(date)] $D R3 pigz done"
else
  echo "[$(date)] $D R3.fastq already gone (already compressed?)"
fi
# integrity check: all three .gz valid
for r in R1 R2 R3; do "$PIGZ" -t "${D}_S1_L001_${r}_001.fastq.gz" && echo "$D $r OK"; done
ls -la
