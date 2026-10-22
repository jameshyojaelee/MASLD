#!/bin/bash
#SBATCH --job-name=V2scr
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=90:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

# V2 target screen: one library per job (pilot). Never submitted as an array.
#
# Usage: sbatch run_screen_one.sbatch <SRR> <COHORT> <RESULT_DIR> [FQ1] [FQ2]

set -euo pipefail

SCRIPTS=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/hypotheses_2026q4/V2_viral_screen

echo "node=$(hostname) job=${SLURM_JOB_ID:-na} started=$(date -u +%FT%TZ)"
bash "$SCRIPTS/screen_library.sh" "$@"
echo "finished=$(date -u +%FT%TZ)"
