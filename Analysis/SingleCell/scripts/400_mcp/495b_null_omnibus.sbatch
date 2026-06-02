#!/bin/bash
#SBATCH --job-name=mcp_null_omnibus
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --mem=100G
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/495b_null_omnibus_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/495b_null_omnibus_%j.err

set -euo pipefail
export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${MASLD_PROJECT_ROOT}/Analysis/SingleCell/scripts/400_mcp"

source "$HOME/.bashrc"
set +u
micromamba activate rnaseq
set -u

echo "[495b.sbatch] host=$(hostname) jobid=${SLURM_JOB_ID:-NA}"
# A10 already delivered by Echo (reviewer_defense/drug_enrichment_null_k16_hvg.tsv); skip.
python 495b_null_omnibus.py --skip a10
echo "[495b.sbatch] DONE."
