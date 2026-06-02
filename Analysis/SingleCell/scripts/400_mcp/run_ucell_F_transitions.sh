#!/bin/bash
#SBATCH --job-name=ucell_Ftrans
#SBATCH --partition=bigmem
#SBATCH --qos=interactive
#SBATCH --mem=500G
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/ucell_Ftrans.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/ucell_Ftrans.err

set -eo pipefail
export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${MASLD_PROJECT_ROOT}/Analysis/SingleCell/scripts/400_mcp"

source "$HOME/.bashrc"
micromamba activate rnaseq

set -u  # defer until after env activation (binutils activate script has unbound vars)

echo "[ucell] host=$(hostname)  cwd=$(pwd)"
echo "[ucell] R=$(which R)"
R --version | head -1

Rscript run_ucell_F_transitions.R

echo "[ucell] DONE."
