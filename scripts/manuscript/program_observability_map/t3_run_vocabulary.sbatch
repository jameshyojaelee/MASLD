#!/bin/bash
#SBATCH --job-name=t3vocab
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/workstreams/BULK-PROGRAM-MAP-v9/logs/t3_vocabulary_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/workstreams/BULK-PROGRAM-MAP-v9/logs/t3_vocabulary_%j.err

set -euo pipefail

RSCRIPT=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
SRC=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/manuscript/program_observability_map

export T3_CORES=16
export T3_RANDOM_DRAWS=${T3_RANDOM_DRAWS:-200}
export T3_CALIBRATION_REPS=${T3_CALIBRATION_REPS:-100}
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

"$RSCRIPT" "$SRC/t3_run_tests.R"
"$RSCRIPT" "$SRC/t3_run_vocabulary.R"
