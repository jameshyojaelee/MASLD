#!/bin/bash
#SBATCH --job-name=ml_target_pred
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=RNA-seq/results/multi_evidence/logs/ml_pred_%j.out
#SBATCH --error=RNA-seq/results/multi_evidence/logs/ml_pred_%j.err

echo "=== ML Target Prediction ==="
echo "Job: $SLURM_JOB_ID | Node: $(hostname) | CPUs: $SLURM_CPUS_PER_TASK"
echo "Started: $(date)"

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

eval "$(micromamba shell hook -s bash)"
micromamba activate rapids_singlecell

export PYTHONUNBUFFERED=1
python -u RNA-seq/45b_ml_target_prediction.py

echo "Finished: $(date)"
