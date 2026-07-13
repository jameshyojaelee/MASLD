#!/usr/bin/env bash
#SBATCH --job-name=saturation
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/102_saturation_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/102_saturation_%j.err
# Submit only after the full array is known: sbatch --dependency=afterok:<ARRAY_JOB_ID> this_file
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export MASLD_PROJECT_ROOT="$ROOT" R_DATATABLE_NUM_PROCS="${SLURM_CPUS_PER_TASK:-8}"
mkdir -p "$ROOT/GWAS/finemapping/results/seqfunc/logs"
micromamba run -n rnaseq Rscript "$ROOT/GWAS/finemapping/src/102_saturation_score_aggregate.R"
