#!/bin/bash
#SBATCH --job-name=causal_chain
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/multi_evidence/logs/causal_chain_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/multi_evidence/logs/causal_chain_%j.err

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "$BASE/RNA-seq"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
set +u
micromamba activate rnaseq
set -u

echo "=== Causal chain (75 -> 217) post-LFC0.5 atlas rebuild ==="
echo "Started: $(date)"

echo ""
echo "--- 75 integrate_causal_overhaul ---"
Rscript 75_integrate_causal_overhaul.R 2>&1

echo ""
echo "--- 217 stratified_causal_atlas ---"
Rscript 217_stratified_causal_atlas.R 2>&1

echo ""
echo "=== DONE at $(date) ==="
