#!/bin/bash
#SBATCH --job-name=masld_donor_context
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=20:00:00
#SBATCH --dependency=afterok:21773588:21773589
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_donor_context_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_donor_context_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 JAX_PLATFORMS=cpu
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PY=/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_local/bin/python
BASE=$PROJ/GWAS/finemapping/results/alphagenome_campaign/week1-20260915
"$PY" "$PROJ/scripts/analysis/alphagenome_campaign/model_donor_tests.py" \
 --out "$BASE/model/donor_fit_invariants_${SLURM_JOB_ID}.json"
"$PY" "$PROJ/scripts/analysis/alphagenome_campaign/model_donor_fit.py" \
 --fixture "$BASE/model/donor_fixture_21773589" \
 --sequence-features "$BASE/model/donor_features_21773588/features.npz" \
 --exclusion-levels chromosome --max-core-hours 79 --cpus 4 --bootstrap 1000 \
 --out "$BASE/model/donor_context_${SLURM_JOB_ID}"
