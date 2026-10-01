#!/bin/bash
#SBATCH --job-name=masld_ag_native_final_review
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_ag_native_final_review_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_ag_native_final_review_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
review_root=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
review_base=$review_root/GWAS/finemapping/results/alphagenome_campaign/week1-20260915
review_python=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
"$review_python" "$review_root/scripts/analysis/alphagenome_campaign/review_native_final.py" \
  --final "$review_base/model/final_21775091" --model "$review_base/model" \
  --sweep "$review_base/model/sweep_21772726" --out "$review_base/review/native-final-$SLURM_JOB_ID"
cd "$review_root"
python3 scripts/manuscript/validate_resource_scope.py
