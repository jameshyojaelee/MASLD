#!/bin/bash
#SBATCH --job-name=convergence_rebuild
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=logs/convergence_rebuild_%j.log

set -euo pipefail
export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "$MASLD_PROJECT_ROOT/RNA-seq"
mkdir -p logs

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

echo "=== Step 1/6: Rebuild atlas (27a) ==="
Rscript 27a_assemble_evidence_atlas.R
echo ""

echo "=== Step 2/6: Re-benchmark presets (27b) ==="
Rscript 27b_benchmark_presets.R
echo ""

echo "=== Step 3/6: Cumulative enrichment (46a) ==="
Rscript 46a_cumulative_enrichment.R
echo ""

echo "=== Step 4/6: Bayesian integration (46b) ==="
Rscript 46b_bayesian_integration.R
echo ""

echo "=== Step 5/6: Network propagation (46c) ==="
python 46c_network_propagation.py
echo ""

echo "=== Step 6/6: Generate Figure 5 ==="
cd "$MASLD_PROJECT_ROOT"
Rscript scripts/figures/fig5_convergence.R
echo ""

echo "=== DONE ==="
ls -la RNA-seq/results/multi_evidence/multi_evidence_atlas.csv
ls -la RNA-seq/results/multi_evidence/bayesian_posterior.csv
ls -la RNA-seq/results/multi_evidence/cumulative_enrichment.csv
ls -la RNA-seq/results/multi_evidence/mutual_information.csv
ls -la figures/fig5/fig5_convergence.pdf
