#!/bin/bash
#SBATCH --job-name=reads_arr
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --array=1-1466%48
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=90:00:00
#SBATCH --output=logs/reads_extraction/arr-%A_%a.out
#SBATCH --error=logs/reads_extraction/arr-%A_%a.err
# One library per task from RUN_DIR/library_manifest.tsv (1,466 libraries with a STAR BAM, largest first).
# NOT SUBMITTED by the agent that wrote it. Usage: sbatch run_extract_array.sbatch RUN_DIR
exec bash /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/hypotheses_2026q4/reads_extraction/extract_task.sh "${1:?run dir}" "$SLURM_ARRAY_TASK_ID"
