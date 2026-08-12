#!/bin/bash
#SBATCH --job-name=masld-fn-map-preflight
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/workstreams/BULK-PROGRAM-MAP-v8/logs/preflight_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/workstreams/BULK-PROGRAM-MAP-v8/logs/preflight_%j.err

set -euo pipefail

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SCRIPT_ROOT=${ROOT}/scripts/manuscript/fibrosis_nas_map
R_BIN=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript

cd "${ROOT}"
for file in "${SCRIPT_ROOT}"/*.R; do
  "${R_BIN}" -e 'parse(file=commandArgs(TRUE)[1])' "${file}" >/dev/null
done
for file in "${SCRIPT_ROOT}"/*.sbatch; do
  bash -n "${file}"
done
python3 - <<'PY'
import ast
from pathlib import Path
for path in Path("scripts/manuscript/fibrosis_nas_map").glob("*.py"):
    ast.parse(path.read_text(), filename=str(path))
print("Python syntax checks passed")
PY
"${R_BIN}" "${SCRIPT_ROOT}/run_tests.R"
