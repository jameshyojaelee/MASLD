#!/bin/bash
#SBATCH --job-name=trail4
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --time=90:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/TRAIL4-20261111T025501Z/logs/trail4_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/TRAIL4-20261111T025501Z/logs/trail4_%j.err
# TRAIL4: one job, no array. Prespecification: ../PRESPEC_TRAIL4.md
set -euo pipefail
RUN=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/TRAIL4-20261111T025501Z
export PATH=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin:$PATH
module load PLINK/1.9
OUT=$RUN/out_${SLURM_JOB_ID}
test ! -e "$OUT"; mkdir -p "$OUT"
echo "host $(hostname) | $(date) | $(plink --version | head -1)"
sha256sum $RUN/PRESPEC_TRAIL4.md $RUN/scripts/* $RUN/downloads/*/*.tsv.gz > $OUT/script_and_input_sha256.txt
cd "$OUT"
Rscript $RUN/scripts/01_trail4_coloc.R "$OUT"
echo "done $(date)"
