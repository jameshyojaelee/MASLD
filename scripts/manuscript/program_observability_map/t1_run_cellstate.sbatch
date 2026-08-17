#!/bin/bash
#SBATCH --job-name=cellstate
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/home/jameslee/t1_work/logs/t1_cellstate_%j.out
#SBATCH --error=/gpfs/commons/home/jameslee/t1_work/logs/t1_cellstate_%j.err

set -euo pipefail

PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
RSCRIPT=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript

export MASLD_T1_CELL_COUNTS=${MASLD_T1_CELL_COUNTS:-/gpfs/commons/home/jameslee/t1_work/run_cell_counts.tsv}
export MASLD_T1_PERMUTATIONS=${MASLD_T1_PERMUTATIONS:-2000}

"$RSCRIPT" "$PROJECT_ROOT/scripts/manuscript/program_observability_map/t1_tests.R"
"$RSCRIPT" "$PROJECT_ROOT/scripts/manuscript/program_observability_map/t1_run_cellstate.R"
