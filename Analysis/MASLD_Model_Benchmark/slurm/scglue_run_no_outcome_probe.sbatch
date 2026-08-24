#!/usr/bin/env bash
#SBATCH --job-name=masld-scglue-no-outcome
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=4
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
BASE="${ROOT}/executions/environments/scvi-tools-1.5.0.post1-portable-21064579"
BASE_ARTIFACTS_SHA256=30814d4944e933d999e8c02c5b521432786b1985cc8794f9436f2a44af5ad1af
GLUE="${ROOT}/executions/environments/scglue-v0.4.1-portable-21064673"
GLUE_ARTIFACTS_SHA256=371632648d40f8838f97605442623c9302886d32f84e0282f1167445ac513081
TORCH="${ROOT}/executions/environments/pytorch-2.6.0-cu124-portable-21064747"
TORCH_ARTIFACTS_SHA256=d4112dbc9574e619edbaec0f2c1503c57e0271fe17517b2e121d711ea21f71a3
PYTHON_SITE="${TORCH}/site-packages:${GLUE}/site-packages:${BASE}/site-packages"
TARGET="${ROOT}/executions/scglue-no-outcome-probe-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

preserve_failed_stage() {
  local failed="${TARGET}.failed"
  if [[ -d "${STAGE}" ]]; then
    if [[ -e "${failed}" ]]; then
      printf 'Cannot preserve failed scGLUE stage: %s\n' "${failed}" >&2
      return
    fi
    mv "${STAGE}" "${failed}"
    printf 'Preserved failed scGLUE stage: %s\n' "${failed}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${TARGET}.failed" ]]; then
  printf 'Refusing to overwrite scGLUE probe: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
printf '%s  %s\n' "${BASE_ARTIFACTS_SHA256}" "${BASE}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${GLUE_ARTIFACTS_SHA256}" "${GLUE}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${TORCH_ARTIFACTS_SHA256}" "${TORCH}/ARTIFACTS.json" \
  | sha256sum --check --strict
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${BASE}" "${GLUE}" "${TORCH}" <<'PY'
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
PYTHONPATH="${PYTHON_SITE}" python "${ROOT}/scripts/scglue_no_outcome_probe.py" \
  --output "${STAGE}/probe" \
  > "${STAGE}/validation/probe.stdout.json"
PYTHONPATH="${PYTHON_SITE}" python -m pip freeze --all \
  > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  "${ROOT}/scripts/scglue_no_outcome_probe.py" \
  "${ROOT}/slurm/scglue_run_no_outcome_probe.sbatch" \
  > "${STAGE}/source.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - \
  "${STAGE}" "${BASE_ARTIFACTS_SHA256}" "${GLUE_ARTIFACTS_SHA256}" \
  "${TORCH_ARTIFACTS_SHA256}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "probe/receipt.json").read_text())
models = receipt["models"]
if (
    receipt["status"] != "pass"
    or receipt["scglue"] != "0.4.1"
    or receipt["torch"] != "2.6.0+cu124"
    or receipt["torch_cuda_build"] != "12.4"
    or "L40S" not in receipt["gpu"]
    or receipt["pairing_topology"] != "same_nucleus"
    or not receipt["ordered_same_nucleus_ids_verified"]
    or receipt["gse244832_false_pairing"]
    or receipt["gse281367_false_pairing"]
    or receipt["project_data_read"]
    or receipt["outcomes_read"]
    or receipt["labels_read"]
    or not receipt["latent_integration_only"]
    or receipt["rna_to_atac_prediction_executed"]
    or receipt["sealed_rna_conditioned_atac_eligible"]
    or receipt["champion_claim_allowed"]
    or set(models) != {"scglue", "paired_scglue"}
    or any(not model["fit_executed"] for model in models.values())
    or any(not model["state_dict_resume_executed"] for model in models.values())
    or any(model["external_dill_deserialization_executed"] for model in models.values())
    or any(model["decode_data_executed"] for model in models.values())
):
    raise SystemExit("scGLUE no-outcome runtime contract differs")
freeze_tree(
    root,
    {
        "artifact_class": "scglue_no_outcome_runtime_probe",
        "model_roster": ["scglue", "paired_scglue"],
        "base_runtime_artifacts_sha256": sys.argv[2],
        "scglue_runtime_artifacts_sha256": sys.argv[3],
        "torch_overlay_artifacts_sha256": sys.argv[4],
        "project_data_read": False,
        "outcomes_read": False,
        "latent_integration_only": True,
        "champion_claim_allowed": False,
        "status": "passed",
    },
)
PY

mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'scGLUE no-outcome probe: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
