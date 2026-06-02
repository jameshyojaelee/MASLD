#!/bin/bash
#SBATCH --job-name=chen_downsample
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --array=1-10%5
#SBATCH --output=logs/chen_downsample_%A_%a.out
#SBATCH --error=logs/chen_downsample_%A_%a.err

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

export DOWNSAMPLE_SIZE=118
export DOWNSAMPLE_ITER=$SLURM_ARRAY_TASK_ID

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts

echo "=== Chen downsampling iteration $SLURM_ARRAY_TASK_ID ==="
echo "Target N: $DOWNSAMPLE_SIZE"
echo "Start: $(date)"

Rscript analysis/integration/scripts/15_chen_downsampling.R

echo "End: $(date)"
