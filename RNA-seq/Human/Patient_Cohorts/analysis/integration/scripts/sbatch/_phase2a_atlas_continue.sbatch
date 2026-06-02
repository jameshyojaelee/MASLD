#!/bin/bash
#SBATCH --job-name=p2a_atlas2
#SBATCH --partition=test
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/multi_evidence/logs/p2a_atlas2_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/multi_evidence/logs/p2a_atlas2_%j.err

# Phase 2A continuation — atlas chain from 45a onward.
# 27a already completed in 15719617 (atlas rebuilt 33,943 genes x 303 columns).
# 35 ATAC (snapatac2 env) skipped — known env conflict; S5 columns inherit prior values.
# 40b sc-TWAS skipped — script no longer exists in repo (only 40_drug_target_validation.R).

set -eo pipefail

BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$BASE/RNA-seq"

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

echo "============================================================"
echo "  Phase 2A continuation (45a → 75 → 57 → 217 → 27b)"
echo "  Started: $(date)"
echo "============================================================"

echo
echo "--- 45a: spatial + scRNA + LIANA + sources_active ---"
Rscript 45a_integrate_all_sources.R

echo
echo "--- 75: causal overhaul ---"
Rscript 75_integrate_causal_overhaul.R

echo
echo "--- 57: ncRNA atlas integration ---"
Rscript 57_ncrna_atlas_integration.R

echo
echo "--- 217: stratified causal atlas integration ---"
Rscript 217_stratified_causal_atlas.R

echo
echo "--- 27b: benchmark presets ---"
Rscript 27b_benchmark_presets.R

echo
echo "============================================================"
echo "  Phase 2A complete: $(date)"
echo "  Atlas: $BASE/RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"
echo "============================================================"
