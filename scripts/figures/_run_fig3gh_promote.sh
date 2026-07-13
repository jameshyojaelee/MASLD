#!/bin/bash
#SBATCH --job-name=fig3gh
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/fig3gh_%j.log
BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$BASE"
mkdir -p scripts/figures/logs

# Activate env BEFORE enabling strict mode: the conda activate.d binutils hook
# references ADDR2LINE unguarded, which trips `set -u` (known env quirk).
source ~/.bashrc 2>/dev/null || true
eval "$(micromamba shell hook --shell bash)" 2>/dev/null || true
micromamba activate rnaseq
set -euo pipefail

echo "=== [1/2] ccc_v3_panels.R (regenerate donor-collapsed ccc_trajectories_data.csv) ==="
Rscript scripts/figures/ccc_v3_panels.R

echo ""
echo "=== [2/2] singlecell_module_heatmap.R (Fig 3G donor-level + Fig 3H from donor-collapsed data) ==="
Rscript scripts/figures/singlecell_module_heatmap.R

echo ""
echo "=== DONE ==="
