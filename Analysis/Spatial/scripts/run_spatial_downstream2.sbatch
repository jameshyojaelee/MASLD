#!/bin/bash
#SBATCH --job-name=spatial_downstream2
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/spatial_downstream2_%j.out
#SBATCH --error=logs/spatial_downstream2_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook -s bash)"
micromamba activate spatial

run_script() {
  local name=$1; local cmd=$2
  echo "=== [$(date +%H:%M:%S)] START $name ==="
  if eval "$cmd" > "logs/${name}_${SLURM_JOB_ID}.out" 2>&1; then
    echo "=== [$(date +%H:%M:%S)] DONE  $name ==="
    return 0
  else
    echo "ERROR: $name failed"; return 1
  fi
}

echo "=============================="
echo "Wave 1: 17a (unified integration)"
echo "=============================="
run_script "17a_unified_spatial_integration" "python scripts/17a_unified_spatial_integration.py"

echo "=============================="
echo "Wave 2: 09b + 09c (parallel)"
echo "=============================="
run_script "09b_spatial_gwas_overlay"   "python scripts/09b_spatial_gwas_overlay.py"   &
run_script "09c_niche_drug_targets"     "python scripts/09c_niche_drug_targets.py"      &
wait

echo "=============================="
echo "Wave 3: 06 (integration)"
echo "=============================="
run_script "06_integration" "python scripts/06_integration.py"

echo "=============================="
echo "Wave 4: 07 + 17b (parallel)"
echo "=============================="
run_script "07_spatial_figures"     "python scripts/07_spatial_figures.py"      &
run_script "17b_spatial_main_figure" "python scripts/17b_spatial_main_figure.py" &
wait

echo "=== Spatial downstream2 complete: $(date) ==="
