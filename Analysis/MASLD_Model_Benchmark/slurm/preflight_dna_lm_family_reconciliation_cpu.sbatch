#!/usr/bin/env bash
#SBATCH --job-name=model-check-216
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1
export WANDB_MODE=offline
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python
CONFIG="${ROOT}/config/dna_lm_family_reconciliation.json"
TARGET="${ROOT}/executions/model-check-216-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite family reconciliation: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
export PYTHONPATH="${ROOT}:${ROOT}/src"

"${PYTHON}" -m unittest -v \
  tests.unit.test_preflight_dna_lm_family_reconciliation \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"

"${PYTHON}" "${ROOT}/scripts/preflight_dna_lm_family_reconciliation.py" \
  --root "${ROOT}" \
  --config "${CONFIG}" \
  --output "${STAGE}/reconciliation" \
  > "${STAGE}/validation/preflight.stdout.json" \
  2> "${STAGE}/validation/preflight.stderr.txt"

"${PYTHON}" -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
"${PYTHON}" --version > "${STAGE}/validation/python.version.txt" 2>&1
sha256sum \
  "${CONFIG}" \
  "${ROOT}/scripts/preflight_dna_lm_family_reconciliation.py" \
  "${ROOT}/tests/unit/test_preflight_dna_lm_family_reconciliation.py" \
  "${ROOT}/slurm/preflight_dna_lm_family_reconciliation_cpu.sbatch" \
  > "${STAGE}/source.sha256"

"${PYTHON}" - "${STAGE}" <<'PY'
import json
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree


stage = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads(
    (stage / "reconciliation/receipt.json").read_text(encoding="utf-8")
)
audit = receipt["evaluation_audit"]
if (
    receipt["status"]
    != "pass_reconciled_development_smoke_no_model_promoted"
    or len(receipt["authorities"]) != 33
    or receipt["model_census_rows"] != 4
    or receipt["promotion_summary_rows"] != 4
    or receipt["representation_context_matched"]
    or not receipt["downstream_head_contract_matched"]
    or not receipt["development_summary_metrics_read"]
    or receipt["raw_reporter_outcomes_parsed"]
    or receipt["sealed_assets_read"]
    or receipt["supervised_eQTL_or_ieQTL_head_fit"]
    or receipt["external_evaluation"]
    or receipt["champion_eligible"]
    or receipt["universal_claim_supported"]
    or receipt["clinical_claim_supported"]
    or audit["best_model_by_mean_spearman"] != "nucleotide_transformer"
    or not audit["best_model_is_restricted"]
    or not audit["all_model_gains_over_allele_identity_ridge_are_negative"]
    or audit["three_seed_screening_recommended"]
    or audit["family_promotion"] != "none"
):
    raise SystemExit("DNA language family reconciliation receipt differs")

freeze_tree(
    stage,
    {
        "artifact_class": "dna_language_family_reconciliation",
        "dataset_id": "gse281364",
        "model_count": 4,
        "development_summary_metrics_read": True,
        "raw_reporter_outcomes_parsed": False,
        "sealed_assets_read": False,
        "external_evaluation": False,
        "champion_eligible": False,
        "family_promotion": "none",
        "status": "pass_reconciled_development_smoke_no_model_promoted",
    },
)
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Family reconciliation: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
