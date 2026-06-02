#!/bin/bash
#SBATCH --job-name=figS_singlecell
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=scripts/figures/logs/fig5_%j.log

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${BASE}"
mkdir -p scripts/figures/logs figures

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

# Ensure reshape2 and ggrastr are available
Rscript -e '
for (pkg in c("ggrastr", "reshape2")) {
  if (!requireNamespace(pkg, quietly=TRUE))
    install.packages(pkg, repos="https://cloud.r-project.org", quiet=TRUE)
}
' 2>/dev/null || true

export MASLD_PROJECT_ROOT="${BASE}"

echo "=== Fig 5: Single-Cell ==="
echo "Start: $(date)"
echo "Node: $(hostname)"

Rscript scripts/figures/figS_singlecell.R

echo "=== Done: $(date) ==="
ls -lh figures/figS_singlecell.pdf 2>/dev/null
