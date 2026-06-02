#!/bin/bash
#SBATCH --job-name=mcp_410_prep
#SBATCH --partition=bigmem
#SBATCH --qos=nslab
#SBATCH --mem=400G
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/410_prep_atlas_cnmf.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/410_prep_atlas_cnmf.err

set -euo pipefail
export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${MASLD_PROJECT_ROOT}/Analysis/SingleCell/scripts/400_mcp"

source "$HOME/.bashrc"
micromamba activate rapids_singlecell

echo "[410.sbatch] host=$(hostname)  cwd=$(pwd)"
echo "[410.sbatch] python=$(which python)"
python --version

python 410_prep_atlas_cnmf.py

echo "[410.sbatch] Now running 412 pseudobulk..."
python 412_prep_pseudobulk_dialogue.py

echo "[410.sbatch] Now running 413 normalize_pseudobulk..."
micromamba deactivate
micromamba activate celltype_bio
Rscript 413_normalize_pseudobulk.R

echo "[410.sbatch] DONE."
