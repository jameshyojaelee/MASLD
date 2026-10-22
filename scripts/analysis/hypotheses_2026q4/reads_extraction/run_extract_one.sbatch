#!/bin/bash
#SBATCH --job-name=reads_one
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=90:00:00
#SBATCH --output=logs/reads_extraction/one-%j.out
#SBATCH --error=logs/reads_extraction/one-%j.err
# Single-library run for testing. Usage: sbatch run_extract_one.sbatch RUN_DIR MANIFEST_INDEX
exec bash /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/hypotheses_2026q4/reads_extraction/extract_task.sh "${1:?run dir}" "${2:?index}"
