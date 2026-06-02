#!/bin/bash
#SBATCH --job-name=pb_de_coarse
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/pb_de_coarse_%j.log

set -eo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
export MASLD_PROJECT_ROOT="${BASE}"

cd "${BASE}/Analysis/SingleCell/scripts"
mkdir -p logs

echo "=== Pseudobulk DE Coarse-Stage Contrasts ==="
echo "Start: $(date)"
echo "Node: $(hostname)"
echo "Job: ${SLURM_JOB_ID:-local}"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

set -u

Rscript "${BASE}/Analysis/SingleCell/scripts/run_pseudobulk_de_coarse_stages.R"

echo "=== Done: $(date) ==="
