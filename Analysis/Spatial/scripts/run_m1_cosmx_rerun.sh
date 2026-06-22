#!/bin/bash
#SBATCH --job-name=cosmx
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/run_m1_cosmx_rerun_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/run_m1_cosmx_rerun_%j.err

# M1 (mega-review A6/A7 remediation): re-run the slide-level CosMx DE (42b, already
# edited to emit all-NaN pval_adj under the 3-vs-1 slide design) and then re-run the
# spatial→atlas integration (06) so the atlas
# spatial_govaere2026_cosmx_<ct>_mash_padj columns become all-NaN instead of the
# stale 779 hep / 223 kc pseudoreplicated per-cell "sig" values.
set -euo pipefail

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
MM=/gpfs/commons/home/jameslee/.local/bin/micromamba
cd "$ROOT"

echo "[$(date)] host=$(hostname) job=${SLURM_JOB_ID:-NA}"
echo "[$(date)] STEP 1/2: 42b_govaere2026_cosmx_de.py (slide-level, NaN padj)"
"$MM" run -n spatial python Analysis/Spatial/scripts/42b_govaere2026_cosmx_de.py

echo "[$(date)] STEP 2/2: 06_integration.py (re-map slide-level cols into atlas)"
"$MM" run -n spatial python Analysis/Spatial/scripts/06_integration.py

echo "[$(date)] DONE"
