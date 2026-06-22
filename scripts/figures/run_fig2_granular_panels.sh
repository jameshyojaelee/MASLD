#!/bin/bash
#SBATCH --job-name=limma
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/fig2_granular_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/fig2_granular_%j.err
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=4:00:00

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
FIGS="${BASE}/scripts/figures"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

# strict mode AFTER activation (conda activate.d scripts reference unbound vars)
set -uo pipefail

export MASLD_PROJECT_ROOT="${BASE}"
mkdir -p "${FIGS}/logs"

echo "=== Re-render fig2 granular staging panels: $(date) ==="

PANELS=(
  nas_stage_degs.R
  nas_stage_upset.R
  fib_stage_degs.R
  fib_stage_upset.R
  cascade_degs.R
)

rc=0
for p in "${PANELS[@]}"; do
  echo ""
  echo "--- ${p} ---"
  if Rscript "${FIGS}/${p}" 2>&1; then
    echo "[ok] ${p}"
  else
    echo "[FAIL] ${p}"
    rc=1
  fi
done

echo ""
echo "=== Done (rc=${rc}): $(date) ==="
exit $rc
