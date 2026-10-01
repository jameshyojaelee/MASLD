#!/bin/bash
#SBATCH --job-name=masld_native_residual
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_native_residual_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_native_residual_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
BASE=$PROJ/GWAS/finemapping/results/alphagenome_campaign/week1-20260915/model
"$PY" "$PROJ/scripts/analysis/alphagenome_campaign/model_native_residual_tests.py" \
  --out "$BASE/native_residual_invariants_${SLURM_JOB_ID}.json"
"$PY" "$PROJ/scripts/analysis/alphagenome_campaign/model_native_residual.py" fit \
  --features "$BASE/frozen" --native "$BASE/native/native_1048576.tsv" \
  --labels "$BASE/native/currin_variants.tsv" \
  --resamples 2000 --cpus 4 --max-seconds 6900 \
  --out "$BASE/native_residual_${SLURM_JOB_ID}"
