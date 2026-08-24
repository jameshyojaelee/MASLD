#!/usr/bin/env bash
#SBATCH --job-name=masld-multivi-smoke-prepare
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
SOURCE_ROOT="${ROOT}/executions/gse296875-rna-atac-smoke-20260824T015422Z-21039865/subset"
SOURCE_ARTIFACTS_SHA256=db9e1cb847ec83a61ec79f7a206485db11e323f6b1b54fa4e10bbe4f9728aa3a
SOURCE="${SOURCE_ROOT}/gse296875_rna_atac_smoke_1000.h5"
SOURCE_SHA256=a1c05a432cec6894acae20a129ecceb153f191d9911843c86a90e13cd7f8ce14
RUNTIME="${ROOT}/executions/environments/scvi-tools-1.5.0.post1-portable-21064579"
RUNTIME_ARTIFACTS_SHA256=30814d4944e933d999e8c02c5b521432786b1985cc8794f9436f2a44af5ad1af
SITE_PACKAGES="${RUNTIME}/site-packages"
TARGET="${ROOT}/executions/multivi-smoke-prepared-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

preserve_failed_stage() {
  local failed="${TARGET}.failed"
  if [[ -d "${STAGE}" ]]; then
    if [[ -e "${failed}" ]]; then return; fi
    mv "${STAGE}" "${failed}"
    printf 'Preserved failed MultiVI preparation: %s\n' "${failed}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${TARGET}.failed" ]]; then
  printf 'Refusing to overwrite MultiVI preparation: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
printf '%s  %s\n' "${SOURCE_ARTIFACTS_SHA256}" "${SOURCE_ROOT}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${SOURCE_SHA256}" "${SOURCE}" | sha256sum --check --strict
printf '%s  %s\n' "${RUNTIME_ARTIFACTS_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${SOURCE_ROOT}" "${RUNTIME}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for path in sys.argv[1:]:
    verify_frozen_tree(Path(path))
PY
module purge
module load python/3.12.3-GCCcore-13.3.0
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
PYTHONPATH="${SITE_PACKAGES}" python "${ROOT}/scripts/multivi_prepare_smoke.py" \
  --source "${SOURCE}" \
  --output "${STAGE}/prepared" \
  > "${STAGE}/validation/prepare.stdout.json"
sha256sum \
  "${ROOT}/scripts/multivi_prepare_smoke.py" \
  "${ROOT}/slurm/multivi_prepare_smoke.sbatch" \
  > "${STAGE}/source.sha256"
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "prepared/prepare_receipt.json").read_text())
if (
    receipt["status"] != "pass"
    or receipt["pairing_topology"] != "same_nucleus"
    or receipt["outer_unit"] != "donor"
    or set(receipt["folds"]) != set(map(str, range(5)))
    or not receipt["held_atac_read_in_prepare"]
    or receipt["held_atac_exported_from_prepare"]
    or receipt["fit_predict_source_hdf5_required"]
    or receipt["gse244832_false_pairing"]
    or receipt["gse281367_false_pairing"]
    or receipt["outcomes_read"]
    or receipt["champion_claim_allowed"]
):
    raise SystemExit("MultiVI preparation contract differs")
freeze_tree(root, {
    "artifact_class": "multivi_smoke_prepared_folds",
    "model_id": "multivi",
    "dataset_id": "gse296875",
    "view_id": "gse296875_rna_atac_smoke_1000_v1",
    "pairing_topology": "same_nucleus",
    "source_artifacts_sha256": "db9e1cb847ec83a61ec79f7a206485db11e323f6b1b54fa4e10bbe4f9728aa3a",
    "runtime_artifacts_sha256": "30814d4944e933d999e8c02c5b521432786b1985cc8794f9436f2a44af5ad1af",
    "held_atac_exported": False,
    "outcomes_read": False,
    "champion_claim_allowed": False,
    "status": "passed",
})
PY
mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'MultiVI smoke prepared folds: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
