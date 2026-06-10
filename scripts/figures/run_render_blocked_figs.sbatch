#!/bin/bash
#SBATCH --job-name=Rscript
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/render_blocked_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/render_blocked_%j.err

set -eo pipefail

BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p "${BASE}/scripts/figures/logs"
cd "${BASE}/scripts/figures"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
set +u
micromamba activate rnaseq
set -u

echo "=== render start: $(date) on $(hostname) ==="

echo ""
echo "############ figS_network_stats.R ############"
Rscript figS_network_stats.R || echo "RC=$? figS_network_stats.R returned nonzero"

echo ""
echo "############ figS_sex_panels.R ############"
Rscript figS_sex_panels.R || echo "RC=$? figS_sex_panels.R returned nonzero"

echo ""
echo "=== render end: $(date) ==="
