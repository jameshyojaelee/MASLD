#!/bin/bash
#SBATCH --job-name=masld_donor_curve_probe_gpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_donor_curve_probe_gpu_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_donor_curve_probe_gpu_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export XLA_PYTHON_CLIENT_PREALLOCATE=false HF_HUB_OFFLINE=1
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PY=/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_local/bin/python
BASE=$PROJ/GWAS/finemapping/results/alphagenome_campaign/week1-20260915/model
FIX=$BASE/donor_curve_fixtures_21784771/size_128
# Deliberately tiny core-hour reservation: the script measures 100 steps, writes
# throughput_100steps.json, then refuses. This job exists to get that measurement
# on an L40S, not to fit anything.
"$PY" "$PROJ/scripts/analysis/alphagenome_campaign/model_donor_fit.py" \
  --fixture "$FIX" --sequence-features "$FIX/features.npz" --expect-regions 1152 \
  --exclusion-levels chromosome --max-core-hours 1.0 --cpus 8 --bootstrap 1000 \
  --out "$BASE/donor_curve_probe_gpu_${SLURM_JOB_ID}"
