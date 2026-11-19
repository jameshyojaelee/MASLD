#!/bin/bash
#SBATCH --job-name=coloc_master_refit
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Hypotheses/EQTL-DEDUP-PSD-20261118/logs/build_refit_master_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Hypotheses/EQTL-DEDUP-PSD-20261118/logs/build_refit_master_%j.err
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
RUN=$ROOT/Analysis/Hypotheses/EQTL-DEDUP-PSD-20261118
PY=/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python
echo "host $(hostname) | $(date -u)"
$PY -I $RUN/code/build_refit_master_20261118.py
sha256sum $RUN/code/build_refit_master_20261118.py \
  $RUN/results/susie_coloc_all_gwas_refit_20261118.csv \
  $RUN/results/r2_summary_refit_20261118.json \
  $RUN/results/step1_counts_20261118.json \
  > $RUN/results/build_refit_master_20261118.sha256
echo "DONE $(date -u)"
