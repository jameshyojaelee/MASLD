#!/usr/bin/env bash
#SBATCH --job-name=masld-lsi-no-outcome
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME="${ROOT}/executions/environments/scvi-tools-1.5.0.post1-portable-21064579"
RUNTIME_ARTIFACTS_SHA256=30814d4944e933d999e8c02c5b521432786b1985cc8794f9436f2a44af5ad1af
SITE_PACKAGES="${RUNTIME}/site-packages"
TARGET="${ROOT}/executions/lsi-no-outcome-probe-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

preserve_failed_stage() {
  local failed="${TARGET}.failed"
  if [[ -d "${STAGE}" ]]; then
    if [[ -e "${failed}" ]]; then
      printf 'Cannot preserve failed LSI stage: %s\n' "${failed}" >&2
      return
    fi
    mv "${STAGE}" "${failed}"
    printf 'Preserved failed LSI stage: %s\n' "${failed}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${TARGET}.failed" ]]; then
  printf 'Refusing to overwrite LSI probe: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
printf '%s  %s\n' "${RUNTIME_ARTIFACTS_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${RUNTIME}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

verify_frozen_tree(Path(sys.argv[1]))
PY
module purge
module load python/3.12.3-GCCcore-13.3.0
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
PYTHONPATH="${SITE_PACKAGES}" python "${ROOT}/scripts/lsi_no_outcome_probe.py" \
  --output "${STAGE}/probe" \
  > "${STAGE}/validation/probe.stdout.json"
PYTHONPATH="${SITE_PACKAGES}" python -m pip freeze --all \
  > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  "${ROOT}/scripts/lsi_no_outcome_probe.py" \
  "${ROOT}/slurm/lsi_run_no_outcome_probe.sbatch" \
  > "${STAGE}/source.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - \
  "${STAGE}" "${RUNTIME_ARTIFACTS_SHA256}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "probe/receipt.json").read_text())
if (
    receipt["status"] != "pass"
    or receipt["scikit_learn"] != "1.9.0"
    or not receipt["held_transform_reused_training_idf"]
    or not receipt["held_transform_did_not_refit_svd"]
    or not receipt["query_composition_invariant"]
    or not receipt["repeat_transform_bit_identical"]
    or not receipt["observed_all_zero_atac_rejected"]
    or receipt["structurally_missing_atac_mapped_to_zero"]
    or not receipt["observed_atac_required_at_query"]
    or receipt["rna_conditioned_atac_eligible"]
    or receipt["sealed_rna_only_eligible"]
    or receipt["project_data_read"]
    or receipt["outcomes_read"]
    or receipt["champion_claim_allowed"]
):
    raise SystemExit("LSI no-outcome contract differs")
freeze_tree(
    root,
    {
        "artifact_class": "lsi_no_outcome_runtime_probe",
        "runtime_artifacts_sha256": sys.argv[2],
        "project_data_read": False,
        "outcomes_read": False,
        "observed_atac_development_only": True,
        "champion_claim_allowed": False,
        "status": "passed",
    },
)
PY

mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'LSI no-outcome probe: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
