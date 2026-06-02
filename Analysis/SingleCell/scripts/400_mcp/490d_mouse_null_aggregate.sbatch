#!/bin/bash
#SBATCH --job-name=mcp_mouse_null_agg
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/490d_agg_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/400_mcp/logs/490d_agg_%j.err
set +u
source "$HOME/.bashrc"
micromamba activate rapids_singlecell
set -u

export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${MASLD_PROJECT_ROOT}/Analysis/SingleCell/scripts/400_mcp"

python 490d_mouse_null_aggregate.py --k 16 --score-kind score_donor_mean
echo "[490d.agg] DONE."
