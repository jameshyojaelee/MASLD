#!/bin/bash
#SBATCH --job-name=masld_native_lengths
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_native_lengths_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_native_lengths_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export XLA_PYTHON_CLIENT_PREALLOCATE=false HF_HUB_OFFLINE=1
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
length_root=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
length_base=$length_root/GWAS/finemapping/results/alphagenome_campaign/week1-20260915
length_python=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
length_scoring_python=/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_local/bin/python
"$length_python" "$length_root/scripts/analysis/alphagenome_campaign/campaign_budget_tests.py"
"$length_python" "$length_root/scripts/analysis/alphagenome_campaign/model_native_lengths_tests.py"
"$length_python" "$length_root/scripts/analysis/alphagenome_campaign/model_native_lengths.py" \
  --native "$length_base/model/native" \
  --adaptation "$length_base/model/sweep_summary_21773124/matched_predictions.tsv.gz" \
  --out "$length_base/model/native_lengths_$SLURM_JOB_ID" \
  --scoring-python "$length_scoring_python" --max-hours 6
