#!/bin/bash
#SBATCH --job-name=prediction
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/leakagefree_chain_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/leakagefree_chain_%j.err
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00

# M7 (mega-review A7 117-fix effectuation): the leakage-free prediction chain.
# 117 re-derives per-fold S1-vs-S2 divergence ranking -> divergence_genes_per_fold.csv
# (replacing the full-cohort 100 div_ cols staged in feature_matrix_v2.csv, Mar-30).
# Chain regenerates the v2 feature matrix and the LOCO / multi-output predictors so the
# honest (leakage-free) S2-fate AUROC numbers land downstream.
# 190 (GPU torch CoralLayer) is submitted SEPARATELY on the gpu partition, not here.

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPTS="${BASE}/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate spatial

# strict mode AFTER activation (conda activate.d scripts reference unbound vars)
set -euo pipefail

export MASLD_PROJECT_ROOT="${BASE}"

cd "${SCRIPTS}"

echo "=== Leakage-free prediction chain: $(date) ==="
echo "SLURM_JOB_ID: ${SLURM_JOB_ID}"
echo "python: $(which python)"
echo ""

echo "--- 117_bifurcation_divergence.py ---"
python 117_bifurcation_divergence.py

echo "--- 150_prognosis_feature_engineering.py ---"
python 150_prognosis_feature_engineering.py

echo "--- 160_feature_engineering_v2.py ---"
python 160_feature_engineering_v2.py

echo "--- 161_nested_loco_classifier.py ---"
python 161_nested_loco_classifier.py

echo "--- 183_multi_output_predictor.py ---"
python 183_multi_output_predictor.py

echo "--- 183v2_multi_output_ablation.py ---"
python 183v2_multi_output_ablation.py

echo ""
echo "=== Done: $(date) ==="
