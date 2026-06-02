#!/bin/bash
#SBATCH --job-name=proteo_panels
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=logs/proteo_panels_%j.out
#SBATCH --error=logs/proteo_panels_%j.err

set -euo pipefail

echo "=== Proteomics Presentation Panels ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Start: $(date)"

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p scripts/figures/logs

eval "$(micromamba shell hook --shell=bash)"
set +u  # conda activation references unbound vars (ADDR2LINE) — aborts under set -u
micromamba activate rnaseq  # matplotlib/pandas + adjustText (CPU); not in spatial, was rapids_singlecell
set -u
export MPLBACKEND=Agg

python scripts/figures/proteomics_presentation_panels.py

echo "=== DONE: $(date) ==="
