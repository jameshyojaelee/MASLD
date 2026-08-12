#!/bin/bash
#SBATCH --job-name=masld-fn-map-discovery
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/workstreams/BULK-PROGRAM-MAP-v8/logs/discovery_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/workstreams/BULK-PROGRAM-MAP-v8/logs/discovery_%j.err

set -euo pipefail

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SCRIPT_ROOT=${ROOT}/scripts/manuscript/fibrosis_nas_map
R_BIN=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript

cd "${ROOT}"
python3 "${SCRIPT_ROOT}/freeze_contract.py"
"${R_BIN}" "${SCRIPT_ROOT}/partition_inputs.R"
"${R_BIN}" "${SCRIPT_ROOT}/run_tests.R"
"${R_BIN}" "${SCRIPT_ROOT}/run_discovery.R"
