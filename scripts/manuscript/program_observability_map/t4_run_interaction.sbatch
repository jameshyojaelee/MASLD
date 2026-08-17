#!/bin/bash
#SBATCH --job-name=t4interaction
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=20:00:00
#SBATCH --output=/gpfs/commons/home/jameslee/t4_run_logs/t4_interaction_%j.log
set -euo pipefail
R=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
# One BLAS thread per worker; the parallelism is over permutation replicates.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export T4_CORES=16
# 140,000 empirical permutations put the per-edge resolution floor at 7.14e-6,
# below the 8.07e-6 BH threshold a 6,195-pair family needs.
export T4_REFERENCE_PERMUTATIONS=4000
export T4_CALIBRATION_PERMUTATIONS=2000
export T4_EMPIRICAL_PERMUTATIONS=140000
export T4_GROUP_PERMUTATIONS=1000

echo "=== synthetic unit tests ==="
$R scripts/manuscript/program_observability_map/t4_run_tests.R
echo "=== T4 interaction run ==="
$R scripts/manuscript/program_observability_map/t4_run_interaction.R
