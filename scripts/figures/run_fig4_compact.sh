#!/usr/bin/env bash
#SBATCH --job-name=fig4_compact
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/fig4_compact_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/fig4_compact_%j.err
set -euo pipefail

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${MASLD_PROJECT_ROOT}"

# Install ggrastr if missing
Rscript -e 'if (!requireNamespace("ggrastr", quietly=TRUE)) install.packages("ggrastr", repos="https://cloud.r-project.org", quiet=TRUE)' 2>/dev/null || true

echo "=== Fig 4 Pharma panels ==="
Rscript scripts/figures/fig4_pharma_panels.R
echo ""

echo "=== Fig 4 Proteomics panels ==="
Rscript scripts/figures/fig4_proteomics_panels.R
echo ""

echo "=== Fig 4 Spatial panels ==="
Rscript scripts/figures/fig4_spatial_panels.R
echo ""

echo "=== Fig 4 Compact (assembled) ==="
Rscript scripts/figures/fig4_compact.R
echo ""

echo "=== Done ==="
ls -lh figures/fig4_*.pdf 2>/dev/null
