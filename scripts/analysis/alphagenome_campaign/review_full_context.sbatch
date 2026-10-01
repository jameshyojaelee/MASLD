#!/bin/bash
#SBATCH --job-name=masld_review_full_context
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --dependency=afterok:21773128:21773115
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_review_full_context_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_review_full_context_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
review_root=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
review_base=$review_root/GWAS/finemapping/results/alphagenome_campaign/week1-20260915
review_python=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
"$review_python" "$review_root/scripts/analysis/alphagenome_campaign/review_specialist_identity.py" \
 --out "$review_base/review/specialist-identity-$SLURM_JOB_ID" --verify-production
"$review_python" "$review_root/scripts/analysis/alphagenome_campaign/review_final_sensitivity.py" \
 --out "$review_base/review/final-sensitivity-$SLURM_JOB_ID"
"$review_python" "$review_root/scripts/analysis/alphagenome_campaign/review_full_context.py" \
 --lx2 "$review_base/model/lx2_full_21773115" \
 --figures "$review_root/figures/candidates/regulatory-week1-20260915-v4" \
 --donor "$review_base/review/b2-21771990" \
 --sweep "$review_base/model/sweep_21772726" \
 --sweep-summary "$review_base/model/sweep_summary_21773124" \
 --transfer "$review_base/model/transfer_21773114" \
 --transfer "$review_base/model/transfer_r2_21773128" \
 --features "$review_base/model/reporter_targets_21772928" \
 --out "$review_base/review/full-context-$SLURM_JOB_ID"
