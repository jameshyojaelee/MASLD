#!/bin/bash
#SBATCH --job-name=ggplot2
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=scripts/figures/logs/rerender_conserved_%j.out
#SBATCH --error=scripts/figures/logs/rerender_conserved_%j.err

set -o pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT_DIR="${BASE}/scripts/figures"
mkdir -p "${SCRIPT_DIR}/logs"

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

export MASLD_PROJECT_ROOT="${BASE}"

# Scripts source load_figure_data.R / publication_theme.R relatively; run from scripts/figures/
cd "${SCRIPT_DIR}"

echo "=== [$(date '+%F %T')] figS_conserved_core.R ==="
Rscript figS_conserved_core.R || echo "WARN: figS_conserved_core.R exited non-zero"

echo "=== [$(date '+%F %T')] figS_published_panel_benchmark.R ==="
Rscript figS_published_panel_benchmark.R || echo "WARN: figS_published_panel_benchmark.R exited non-zero"

echo "=== [$(date '+%F %T')] DONE ==="
