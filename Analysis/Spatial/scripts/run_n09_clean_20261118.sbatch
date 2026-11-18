#!/usr/bin/env bash
#SBATCH --job-name=n09_geomx_clean
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=00:30:00
#SBATCH --output=Analysis/Spatial/logs/n09_geomx_clean_%j.out
#SBATCH --error=Analysis/Spatial/logs/n09_geomx_clean_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PY_BIN="/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python"
R_BIN="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"

GUARD_LOG="Analysis/Spatial/logs/n09_guard_${SLURM_JOB_ID}.log"
set +e
"$R_BIN" --vanilla Analysis/Spatial/scripts/40_govaere2026_geomx_load.R > "$GUARD_LOG" 2>&1
rc=$?
set -e
echo "[N09] guard rc=$rc"
cat "$GUARD_LOG"
rm -f "$GUARD_LOG"
if [ "$rc" -eq 0 ]; then
  echo "ERROR: 40_govaere2026_geomx_load.R did not stop!" >&2
  exit 1
fi
echo "[N09] Running 44_spatial_consensus.py & 06_integration.py..."
"$PY_BIN" Analysis/Spatial/scripts/44_spatial_consensus.py
"$PY_BIN" Analysis/Spatial/scripts/06_integration.py
