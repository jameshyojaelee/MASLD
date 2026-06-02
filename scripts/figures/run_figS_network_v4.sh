#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --time=48:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --job-name=figS_network_v4
#SBATCH --output=logs/figS_network_v4_%j.out
#SBATCH --error=logs/figS_network_v4_%j.err
# ============================================================================
# Submit the consolidated v4 supplementary panels (figS_network_v4.R).
# ============================================================================
set -eo pipefail

BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$BASE"
mkdir -p "$BASE/scripts/figures/logs"

# shellcheck disable=SC1091
source "$HOME/.bashrc" || true
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

export MASLD_PROJECT_ROOT="$BASE"

echo "[$(date)] Running figS_network_v4.R (all panels)"
Rscript "$BASE/scripts/figures/figS_network_v4.R"

echo "[$(date)] Done."
