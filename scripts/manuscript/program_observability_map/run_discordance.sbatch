#!/bin/bash
# bigmem variant of run_discovery.sbatch.
#
# This job needs about 96G. The 500G ask is not a real memory requirement: on
# this cluster a bigmem request below 500G is silently redirected to cpu, and
# landing on bigmem also requires --qos=interactive. Both lines below are the
# cost of using the partition at all, not an estimate of the job's footprint.
# --qos=interactive caps the user at 4 concurrent jobs, so this consumes one of
# those slots.
#SBATCH --job-name=v9discord
#SBATCH --partition=bigmem
#SBATCH --qos=interactive
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=500G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/workstreams/BULK-PROGRAM-MAP-v9/logs/discordance_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/workstreams/BULK-PROGRAM-MAP-v9/logs/discordance_%j.err

set -euo pipefail
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
# 50,000 permutations, not a round number chosen for comfort. Per-donor calls
# need an empirical p floor of 1/(B+1) below the BH threshold 0.05/469 = 1.07e-4,
# so B must exceed ~9,400. At B = 1,000 the floor was 9.99e-4 and no donor could
# be called at any effect size, which made the earlier "0 discordant" arithmetic
# rather than biology. 50,000 gives a floor of 2.0e-5, comfortably clear.
export MASLD_DISCORDANCE_PERMUTATIONS=${MASLD_DISCORDANCE_PERMUTATIONS:-50000}

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SCRIPT_ROOT=${ROOT}/scripts/manuscript/program_observability_map
R_BIN=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript

cd "${ROOT}"
"${R_BIN}" "${SCRIPT_ROOT}/run_discordance.R"
