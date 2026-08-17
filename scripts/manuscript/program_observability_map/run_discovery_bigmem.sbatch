#!/bin/bash
# bigmem variant of run_discovery.sbatch.
#
# This job needs about 96G. The 500G ask is not a real memory requirement: on
# this cluster a bigmem request below 500G is silently redirected to cpu, and
# landing on bigmem also requires --qos=interactive. Both lines below are the
# cost of using the partition at all, not an estimate of the job's footprint.
# --qos=interactive caps the user at 4 concurrent jobs, so this consumes one of
# those slots.
#SBATCH --job-name=v9discovery
#SBATCH --partition=bigmem
#SBATCH --qos=interactive
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=500G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/workstreams/BULK-PROGRAM-MAP-v9/logs/discovery_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/workstreams/BULK-PROGRAM-MAP-v9/logs/discovery_%j.err

set -euo pipefail
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SCRIPT_ROOT=${ROOT}/scripts/manuscript/program_observability_map
R_BIN=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript

cd "${ROOT}"
"${R_BIN}" "${SCRIPT_ROOT}/run_discovery.R"
