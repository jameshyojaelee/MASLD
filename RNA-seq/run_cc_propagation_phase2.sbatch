#!/bin/bash
#SBATCH --job-name=ccprop
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/cc_propagation_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/cc_propagation_%j.out

# Phase 2 of the Conserved-Core propagation (after the 27a->75->217 atlas rebuild that
# refreshed is_conserved to 1,355). Re-runs the discrete CC-dependent enrichment scripts
# on the fresh atlas. Each step is NON-FATAL so one failure cannot block the others.
set -uo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
set +u; micromamba activate rnaseq; set -u
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$MASLD_PROJECT_ROOT"
echo "host=$(hostname) start=$(date)"

echo "=== [1] evidence_independence_metrics (triple-convergence DEG+COLOC+CC) ==="
Rscript scripts/evidence_independence_metrics.R 2>&1 || echo "WARN: evidence_independence_metrics rc=$?"

echo "=== [2] 46a_cumulative_enrichment (V3 Conserved validation + cumulative ORs) ==="
Rscript RNA-seq/46a_cumulative_enrichment.R 2>&1 || echo "WARN: 46a rc=$?"

echo "=== [3] H6_inv_cc_positive_control (drug/multi-evidence/positive-control ORs, CC protein rho selection) ==="
Rscript RNA-seq/audit/H6_inv_cc_positive_control.R 2>&1 || echo "WARN: H6 rc=$?"

echo "=== [4] mrna_protein_concordance (Conserved_Core protein rho) ==="
Rscript Analysis/Proteomics/scripts/mrna_protein_concordance.R 2>&1 || echo "WARN: proteomics rc=$?"

echo "=== DONE $(date) ==="
