#!/bin/bash
#SBATCH --job-name=mcp_composition_ctrl
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --mem=100G
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/463_composition_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/463_composition_%j.err

set -euo pipefail
export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${MASLD_PROJECT_ROOT}/Analysis/SingleCell/scripts/400_mcp"

source "$HOME/.bashrc"
set +u
micromamba activate rnaseq
set -u

echo "[463.sbatch] host=$(hostname) jobid=${SLURM_JOB_ID:-NA}"
python 463_composition_controls.py --n-shuffle 1000
echo "[463.sbatch] DONE."
