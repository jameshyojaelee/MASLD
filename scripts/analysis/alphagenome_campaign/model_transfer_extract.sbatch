#!/bin/bash
#SBATCH --job-name=masld_ag_reporter_targets
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_ag_reporter_targets_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_ag_reporter_targets_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export XLA_PYTHON_CLIENT_PREALLOCATE=false HF_HUB_OFFLINE=1
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PY=/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_local/bin/python
BASE=$PROJ/GWAS/finemapping/results/alphagenome_campaign/week1-20260915
"$PY" "$PROJ/scripts/analysis/alphagenome_campaign/model_transfer_extract.py" \
 --reporters "$BASE/data/endpoints-revised-21772544/reporter_treatment_replicates.tsv.gz" \
 --pairs "$BASE/data/endpoints-revised-21772544/reporter_caqtl_all_nominal_pairs.tsv.gz" \
 --folds "$BASE/data/folds-21772544/reporter_fold_companion.tsv" \
 --out "$BASE/model/reporter_targets_${SLURM_JOB_ID}" --max-hours 0.95
