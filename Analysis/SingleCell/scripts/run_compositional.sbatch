#!/bin/bash
#SBATCH --job-name=clr_lmm_compositional
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/logs/compositional_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/logs/compositional_%j.err

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"

echo "=== CLR + LMM Compositional Analysis ==="
echo "Job ID:   $SLURM_JOB_ID"
echo "Node:     $SLURMD_NODENAME"
echo "Start:    $(date)"

mkdir -p "${BASE}/Analysis/SingleCell/scripts/logs"

# Activate rnaseq conda environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

Rscript "${BASE}/Analysis/SingleCell/scripts/compositional_clr_lmm.R"

echo ""
echo "=== Done: $(date) ==="
echo "Outputs in: ${BASE}/Analysis/SingleCell/results_gpu_v2/compositional/"
ls -lh "${BASE}/Analysis/SingleCell/results_gpu_v2/compositional/"
