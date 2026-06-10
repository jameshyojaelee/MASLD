#!/bin/bash
#SBATCH --job-name=ccprop
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/cc_propagation_b_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/cc_propagation_b_%j.out

# Phase 2b: re-run ONLY the two scripts that failed in Phase 2 on a stale dream_* column
# name (now migrated to bulk_*). H6 + proteomics already succeeded and are not repeated.
set -uo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
set +u; micromamba activate rnaseq; set -u
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$MASLD_PROJECT_ROOT"
echo "host=$(hostname) start=$(date)"

echo "=== [1] evidence_independence_metrics (triple-convergence) ==="
Rscript scripts/evidence_independence_metrics.R 2>&1 || echo "WARN: evidence_independence_metrics rc=$?"

echo "=== [2] 46a_cumulative_enrichment (V3 Conserved validation + cumulative ORs) ==="
Rscript RNA-seq/46a_cumulative_enrichment.R 2>&1 || echo "WARN: 46a rc=$?"

echo "=== DONE $(date) ==="
