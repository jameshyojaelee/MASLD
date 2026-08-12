#!/bin/bash
#SBATCH --job-name=masld-program-map-compact
#SBATCH --partition=bigmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=500G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/workstreams/BULK-PROGRAM-MAP-v8/logs/compact_map_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/workstreams/BULK-PROGRAM-MAP-v8/logs/compact_map_%j.err

set -euo pipefail

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SCRIPT_ROOT=${ROOT}/scripts/manuscript/fibrosis_nas_map
R_BIN=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
WORKSTREAM=${ROOT}/RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/workstreams/BULK-PROGRAM-MAP-v8
OUTPUT=${WORKSTREAM}/figures-impact-v5

cd "${ROOT}"
"${R_BIN}" -e 'parse(file=commandArgs(TRUE)[1])' "${SCRIPT_ROOT}/render_program_map_compact.R" >/dev/null
python3 - "${SCRIPT_ROOT}/validate_program_map_compact.py" <<'PY'
import ast
import pathlib
import sys
path = pathlib.Path(sys.argv[1])
ast.parse(path.read_text(), filename=str(path))
PY
"${R_BIN}" "${SCRIPT_ROOT}/render_program_map_compact.R" "${WORKSTREAM}" "${OUTPUT}"
python3 "${SCRIPT_ROOT}/validate_program_map_compact.py" "${WORKSTREAM}" "${OUTPUT}"
