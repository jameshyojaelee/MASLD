#!/bin/bash
#SBATCH --job-name=ccresid
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/cc_residuals_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/cc_residuals_%j.out

# Final residuals of the Conserved-Core propagation: the two enrichment ORs whose
# producers were not in the Phase-2 batch. Both read the fresh atlas/concordance
# (is_conserved=1,355). Non-fatal between steps.
set -uo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
set +u; micromamba activate rnaseq; set -u
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$MASLD_PROJECT_ROOT"
echo "host=$(hostname) start=$(date)"

echo "=== [1] 43_zonation_classification (PP_vs_ConservedCore OR) ==="
Rscript RNA-seq/43_zonation_classification.R 2>&1 || echo "WARN: 43 rc=$?"

echo "=== [2] 33_audit_sensitivity_analyses Part D (CC pathway/drug/top5 ORs) ==="
Rscript RNA-seq/33_audit_sensitivity_analyses.R 2>&1 || echo "WARN: 33 rc=$?"

echo "=== DONE $(date) ==="
