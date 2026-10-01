#!/bin/bash
#SBATCH --job-name=masld_profile_input_stability
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_profile_input_stability_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_profile_input_stability_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
OUT=$PROJ/GWAS/finemapping/results/alphagenome_campaign/two-models-20260922T203900EDT/profile_input_stability_${SLURM_JOB_ID}
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python \
  "$PROJ/scripts/analysis/masld_two_models/profile_input_stability.py" --out "$OUT"
