#!/bin/bash
#SBATCH --job-name=figS_proteo
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=logs/figS_proteo_%j.out
#SBATCH --error=logs/figS_proteo_%j.err

set -euo pipefail

echo "=== Supplementary Proteomics Validation Panels (v3) ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Start: $(date)"

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p scripts/figures/logs

eval "$(micromamba shell hook --shell=bash)"
set +u  # conda activation references unbound vars (ADDR2LINE) — aborts under set -u
micromamba activate spatial  # matplotlib/pandas only — no GPU/scanpy (was rapids_singlecell)
set -u
export MPLBACKEND=Agg

python scripts/figures/figS_proteomics_validation.py

echo "=== DONE: $(date) ==="
