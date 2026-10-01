#!/bin/bash
#SBATCH --job-name=masld_donor_curve_fixtures
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_donor_curve_fixtures_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_donor_curve_fixtures_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 JAX_PLATFORMS=cpu
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PY=/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_local/bin/python
BASE=$PROJ/GWAS/finemapping/results/alphagenome_campaign/week1-20260915/model
cd "$PROJ/scripts/analysis/alphagenome_campaign"
"$PY" model_donor_curve_fixtures.py \
  --panel "$BASE/donor_curve_panel_21783074" \
  --sequence-features "$BASE/donor_curve_features_21784722/features.npz" \
  --out "$BASE/donor_curve_fixtures_${SLURM_JOB_ID}"
