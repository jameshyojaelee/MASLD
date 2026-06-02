#!/bin/bash
#SBATCH --job-name=352b_fstage3way
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --output=Analysis/SingleCell/results_gpu_v2_phase05/logs/352b_fstage3way_%j.out
#SBATCH --error=Analysis/SingleCell/results_gpu_v2_phase05/logs/352b_fstage3way_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p Analysis/SingleCell/results_gpu_v2_phase05/logs

echo "[$(date)] starting 352b_fstage_on_latent (3 latents)"
for L in X_scVI X_harmony X_scanorama; do
  echo "[$(date)] -> latent=$L"
  micromamba run -n rapids_singlecell python -u \
      Analysis/SingleCell/scripts/v2/352b_fstage_on_latent.py --latent $L
done
echo "[$(date)] DONE"
