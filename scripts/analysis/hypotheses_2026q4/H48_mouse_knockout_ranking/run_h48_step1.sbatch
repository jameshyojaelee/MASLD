#!/bin/bash
#SBATCH --job-name=H48_step1
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --time=90:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --output=%x_%j.log
set -euo pipefail
OUT=$1
R=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
PY=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/.envs/perturbation_gears/bin/python
export PYTHONNOUSERSITE=1
$PY -I $R/scripts/analysis/hypotheses_2026q4/H48_mouse_knockout_ranking/h48_step1_prepare.py \
  --h5ad $R/data/perturbation/datasets/saunders2025/processed/saunders_perturb_hep_confident.h5ad \
  --orth /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/external/orthologs/mouse_human_orthologs_symbols.tsv.gz \
  --sig_sc $R/Analysis/SingleCell/results_gpu_v2/pseudobulk_de_donorcollapsed/Hepatocytes_de.csv \
  --sig_bulk $R/RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv \
  --h5dir $R/data/perturbation/datasets/saunders2025/raw/h5 \
  --out $OUT
$PY -m pip freeze > $OUT/pip_freeze_step1.txt
