#!/bin/bash
#SBATCH --job-name=conv_step5_6
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/conv_step5_6_%j.log

set -euo pipefail
export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "$MASLD_PROJECT_ROOT/RNA-seq"

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

echo "=== Step 5/6: Network propagation (46c) ==="
# Use rapids_singlecell env for working pandas (rnaseq env has ABI mismatch)
micromamba activate rapids_singlecell
PYTHONNOUSERSITE=1 python 46c_network_propagation.py
micromamba activate rnaseq
echo ""

echo "=== Step 6/6: Generate Figure 5 ==="
cd "$MASLD_PROJECT_ROOT"
Rscript scripts/figures/fig5_convergence.R
echo ""

echo "=== DONE ==="
ls -la RNA-seq/results/multi_evidence/network_propagation_scores.csv
ls -la figures/fig5/fig5_convergence.pdf
