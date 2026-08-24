#!/usr/bin/env bash
#SBATCH --job-name=masld-mpra-recon-build
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SOURCE="${ROOT}/scripts/gse281364_reconstruct_oligos.cpp"
VALIDATOR="${ROOT}/scripts/validate_gse281364_reconstruction.py"
TARGET="${ROOT}/executions/tools/gse281364-reconstruction-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
    printf 'Preserved failed MPRA reconstruction build: %s\n' "${FAILED}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite reconstruction tool artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/bin" "${STAGE}/validation"

module load GCC/13.2.0
g++ --version > "${STAGE}/validation/gxx.version.txt"
g++ -std=c++17 -O3 -Wall -Wextra -Wpedantic \
  "${SOURCE}" -lz -o "${STAGE}/bin/gse281364_reconstruct_oligos"
ldd "${STAGE}/bin/gse281364_reconstruct_oligos" \
  > "${STAGE}/validation/linked_libraries.txt"
module purge

"${PROJECT}/.agent/scanpy_env/bin/python" - \
  "${VALIDATOR}" "${STAGE}/validation/python_validation.json" <<'PY'
import ast
import json
from pathlib import Path
import sys

source = Path(sys.argv[1]).read_text(encoding="utf-8")
ast.parse(source, filename=sys.argv[1])
receipt = {
    "validator_ast_parse": True,
    "reference_reader": "internal_uncompressed_fixed_width_fasta_index",
    "validator_executed": False,
    "observed_outcomes_loaded": False,
    "sealed_outcomes_loaded": False,
}
Path(sys.argv[2]).write_text(json.dumps(receipt, sort_keys=True, indent=2) + "\n")
PY
"${PROJECT}/.agent/scanpy_env/bin/python" -m pip freeze \
  > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  "${SOURCE}" "${VALIDATOR}" \
  "${ROOT}/slurm/gse281364_build_reconstruction_tool.sbatch" \
  > "${STAGE}/source.sha256"

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import os
import sys

from masld_bench.artifacts import freeze_tree
from masld_bench.hashing import sha256_file

root = Path(sys.argv[1]).resolve(strict=True)
binary = root / "bin/gse281364_reconstruct_oligos"
probe = json.loads((root / "validation/python_validation.json").read_text())
if (
    not binary.is_file()
    or not os.access(binary, os.X_OK)
    or not probe["validator_ast_parse"]
    or probe["validator_executed"]
    or probe["observed_outcomes_loaded"]
    or probe["sealed_outcomes_loaded"]
):
    raise SystemExit("MPRA reconstruction build receipt differs")
freeze_tree(
    root,
    {
        "artifact_class": "gse281364_reconstruction_tool",
        "language": "C++17",
        "binary_sha256": sha256_file(binary),
        "paired_fastq_passes": 2,
        "observed_oligo_bp": 107,
        "barcode_bp": 20,
        "observed_outcomes_loaded": False,
        "sealed_outcomes_loaded": False,
    },
)
PY
module purge

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'GSE281364 reconstruction tool: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
