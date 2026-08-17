#!/bin/bash
#SBATCH --job-name=t6modifier
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/t6_modifiers_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/t6_modifiers_%j.err

set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
RSCRIPT=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
cd "$ROOT"

echo "=== unit tests: planted interaction must be detected before a null is reported ==="
"$RSCRIPT" scripts/manuscript/program_observability_map/t6_run_tests.R

echo "=== T6 sex, age and cohort effect modification ==="
"$RSCRIPT" scripts/manuscript/program_observability_map/t6_run_modifiers.R
