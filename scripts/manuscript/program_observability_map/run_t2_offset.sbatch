#!/bin/bash
#SBATCH --job-name=t2offset
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/workstreams/BULK-PROGRAM-MAP-v9/logs/t2_offset_%j.out
set -euo pipefail
SCRIPTS=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/manuscript/program_observability_map
R=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
export MASLD_T2_PERMUTATIONS=2000
export OMP_NUM_THREADS=4
"$R" "$SCRIPTS/t2_unit_tests.R"
"$R" "$SCRIPTS/t2_donor_offset_diagnostic.R"
