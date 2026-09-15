#!/usr/bin/env bash
#SBATCH --job-name=model-check-072
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_ARTIFACTS_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
TARGET=${ROOT}/executions/model-check-072-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite cross-cohort authority validation.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
cd "${ROOT}"
printf '%s  %s\n' "${RUNTIME_ARTIFACTS_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
export PYTHONPATH="${ROOT}/src:${ROOT}"
"${PYTHON}" -m unittest -v \
  tests.contract.test_cross_cohort_plan \
  tests.contract.test_gse274114_activation_contract \
  tests.unit.test_activate_gse274114_no_fit \
  tests.unit.test_build_gse105127_assay_native_activation \
  tests.unit.test_gse268273_raw_campaign \
  > "${STAGE}/tests.stdout.txt" 2> "${STAGE}/tests.stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum \
  OVERALL_PLAN.md \
  README.md \
  config/datasets/gse105127_zonated_rna_rrbs.toml \
  config/datasets/gse274114_mash_hbv.toml \
  config/evaluation/cross_cohort_expansion.toml \
  config/evaluation/gse105127_adjacent_section_rna_rrbs_task.toml \
  config/evaluation/gse274114_within_instrument_etiology_task.toml \
  config/evaluation/gse274114_within_instrument_etiology_promotion_gate.json \
  scripts/activate_gse274114_no_fit.py \
  scripts/build_gse105127_assay_native_activation.py \
  scripts/plan_gse268273_raw_campaign.py \
  tests/contract/test_cross_cohort_plan.py \
  tests/contract/test_gse274114_activation_contract.py \
  tests/unit/test_activate_gse274114_no_fit.py \
  tests/unit/test_build_gse105127_assay_native_activation.py \
  tests/unit/test_gse268273_raw_campaign.py \
  slurm/validate_cross_cohort_authorities_cpu.sbatch \
  executions/model-data-080-21082249-gse274114/ARTIFACTS.json \
  executions/gse105127-assay-native-activation-21082444/ARTIFACTS.json \
  executions/gse105127-assay-native-activation-21082444/fixture/ARTIFACTS.json \
  executions/model-data-075-21082574/ARTIFACTS.json \
  executions/model-data-075-21082574/plan/ARTIFACTS.json \
  "${RUNTIME}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
freeze_tree(
    root,
    {
        "artifact_class": "cross_cohort_authority_validation",
        "cohort_families": [
            "gse105127_zonated_rna_rrbs",
            "gse268273_imid_masld",
            "gse274114_mash_hbv",
        ],
        "outcomes_read": False,
        "fit_or_score_performed": False,
        "status": "passed",
    },
)
verify_frozen_tree(root)
PY
mv "${STAGE}" "${TARGET}"
trap - EXIT
sha256sum "${TARGET}/ARTIFACTS.json"
