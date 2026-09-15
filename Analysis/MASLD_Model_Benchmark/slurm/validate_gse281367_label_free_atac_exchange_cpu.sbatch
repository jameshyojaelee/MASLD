#!/usr/bin/env bash
#SBATCH --job-name=model-check-996
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON="${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python"
TARGET="${ROOT}/executions/gse281367-label-free-atac-exchange-validation-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite label-free exchange validation: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"

for pair in \
  "5329086965d81bfab7a616e1fa423754da614a11fbc95c452520aba3a0f42237 config/evaluation/gse281367_label_free_atac_exchange.toml" \
  "7c9e5dcc4de2d43e6785ef84599284be67ce7f52f30e94d492e264ed20c13b5f scripts/build_gse281367_label_free_atac_exchange.py" \
  "3705dfbf632b7585649c1b78c3d26cc8682d07156bcdd1534300ac55cf0f2dd1 tests/unit/test_gse281367_label_free_atac_exchange.py" \
  "30fec3c271235b1348026a9c3ceaabe8dd97792dbf774b5abb4a8d5a5cd3f789 scripts/commit_gse281367_atac_transport_predictions.py" \
  "f442ef9a25993e7ccfaff98025056946900def866fad969856c1a71f9c9890c3 tests/unit/test_commit_gse281367_atac_transport_predictions.py" \
  "079a75b9ff59a62246ea0c1caf8eb60686d15e1262a272885a12789a53c3ba9e scripts/score_gse281367_atac_transport.py" \
  "e34c0a35a14a4a49f10a606925bb95f000da3e87c067f631ceb4b40f4dfb6805 tests/unit/test_score_gse281367_atac_transport.py"; do
  read -r expected relative <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${ROOT}/${relative}" | sha256sum --check --strict
done

PYTHONPATH="${ROOT}/src:${ROOT}" "${PYTHON}" -m unittest -v \
  tests.unit.test_gse281367_label_free_atac_exchange \
  tests.unit.test_commit_gse281367_atac_transport_predictions \
  tests.unit.test_score_gse281367_atac_transport \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
"${PYTHON}" -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
"${PYTHON}" --version > "${STAGE}/validation/python.version.txt" 2>&1
sha256sum \
  config/evaluation/gse281367_label_free_atac_exchange.toml \
  scripts/build_gse281367_label_free_atac_exchange.py \
  tests/unit/test_gse281367_label_free_atac_exchange.py \
  scripts/commit_gse281367_atac_transport_predictions.py \
  tests/unit/test_commit_gse281367_atac_transport_predictions.py \
  scripts/score_gse281367_atac_transport.py \
  tests/unit/test_score_gse281367_atac_transport.py \
  slurm/validate_gse281367_label_free_atac_exchange_cpu.sbatch \
  > "${STAGE}/source.sha256"

PYTHONPATH="${ROOT}/src" "${PYTHON}" - "${STAGE}" "${SLURM_JOB_ID}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
stderr = (root / "validation/tests.stderr.txt").read_text(encoding="utf-8")
if "Ran 15 tests" not in stderr or not stderr.rstrip().endswith("OK"):
    raise SystemExit("label-free exchange tests did not pass exactly")
freeze_tree(
    root,
    {
        "artifact_class": "gse281367_label_free_atac_exchange_validation",
        "slurm_job_id": sys.argv[2],
        "tests": 15,
        "label_free_exchange_contract_valid": True,
        "prediction_commit_registered": True,
        "condition_blind_scorer_registered": True,
        "sequence_and_observed_atac_families_separate": True,
        "condition_values_read": False,
        "biological_outcome_arrays_read": False,
        "metrics_calculated": False,
        "champion_claim_allowed": False,
        "status": "passed",
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Label-free ATAC exchange validation: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
