#!/usr/bin/env bash
#SBATCH --job-name=masld-mpra-sei-mpralegnet-eval
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=2
export MKL_NUM_THREADS=2

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PYTHON="${PROJECT}/.agent/scanpy_env/bin/python"
SEI="${ROOT}/executions/sei-gse281364-smoke-21069658"
SEI_SHA=dd22116a2927712d1c7cafc1a3ca654783efce2173c44f010cd9fc22034c1803
MPRA="${ROOT}/executions/model-work-102-21069877"
MPRA_SHA=a59c82c5093a61025f39d5fb12c61820711f58e001ee9bc7858bd505140ec5f9
OUTCOMES="${ROOT}/executions/gse281364-replicate-outcomes-21066307"
OUTCOMES_SHA=cb42bd4037127a4e4d6a406d6f7b00c706bc93b631cdf93083ce994809e791db
TARGET="${ROOT}/executions/gse281364-sei-mpralegnet-evaluation-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
    printf 'Preserved failed MPRA evaluation stage: %s\n' "${FAILED}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite MPRA evaluation artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"

PYTHONPATH="${ROOT}" "${PYTHON}" -m unittest -v \
  tests.unit.test_evaluate_gse281364_sei_mpralegnet \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"

PYTHONPATH="${ROOT}" "${PYTHON}" \
  "${ROOT}/scripts/evaluate_gse281364_sei_mpralegnet.py" \
  --sei "${SEI}" \
  --sei-artifacts-sha256 "${SEI_SHA}" \
  --mpralegnet "${MPRA}" \
  --mpralegnet-artifacts-sha256 "${MPRA_SHA}" \
  --outcomes "${OUTCOMES}" \
  --outcomes-artifacts-sha256 "${OUTCOMES_SHA}" \
  --output "${STAGE}/evaluation" \
  > "${STAGE}/validation/evaluation.stdout.json" \
  2> "${STAGE}/validation/evaluation.stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  "${ROOT}/scripts/evaluate_gse281364_sei_mpralegnet.py" \
  "${ROOT}/tests/unit/test_evaluate_gse281364_sei_mpralegnet.py" \
  "${ROOT}/slurm/evaluate_gse281364_sei_mpralegnet_cpu.sbatch" \
  > "${STAGE}/source.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import csv
import gzip
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "evaluation/receipt.json").read_text())
with (root / "evaluation/metrics.tsv").open(newline="", encoding="utf-8") as handle:
    metrics = list(csv.DictReader(handle, delimiter="\t"))
with gzip.open(
    root / "evaluation/matched_endpoints.tsv.gz", "rt", newline="", encoding="utf-8"
) as handle:
    matched = list(csv.DictReader(handle, delimiter="\t"))
if (
    receipt["status"] != "pass_descriptive_development_evaluation"
    or receipt["outcome_role"] != "exposed_development_MPRA_only"
    or receipt["sealed_outcomes_read"]
    or receipt["experimental_replicates_per_context"] != 4
    or receipt["replicates_used_as_independent_donors"]
    or receipt["donor_count"] != 0
    or receipt["sei_elements"] != 1033
    or receipt["mpralegnet_elements"] != 4359
    or receipt["locus_sequence_groups"] != 1033
    or receipt["models_forced_to_common_endpoint"]
    or receipt["models_ranked"]
    or not receipt["mpralegnet_secondary_nonchampion"]
    or receipt["champion_eligible"]
    or len(receipt["metric_summary"]) != 6
    or len(metrics) != 36
    or len(matched) != 12850
):
    raise SystemExit("GSE281364 development evaluation receipt differs")
freeze_tree(
    root,
    {
        "artifact_class": "gse281364_sei_mpralegnet_development_evaluation",
        "dataset_id": "gse281364",
        "model_ids": ["sei", "mpralegnet_hepg2_test1_val2"],
        "sei_prediction_artifacts_sha256": "dd22116a2927712d1c7cafc1a3ca654783efce2173c44f010cd9fc22034c1803",
        "mpralegnet_prediction_artifacts_sha256": "a59c82c5093a61025f39d5fb12c61820711f58e001ee9bc7858bd505140ec5f9",
        "outcome_artifacts_sha256": "cb42bd4037127a4e4d6a406d6f7b00c706bc93b631cdf93083ce994809e791db",
        "outcome_role": "exposed_development_MPRA_only",
        "experimental_replicates_per_context": 4,
        "replicates_used_as_independent_donors": False,
        "donor_count": 0,
        "models_forced_to_common_endpoint": False,
        "models_ranked": False,
        "mpralegnet_secondary_nonchampion": True,
        "champion_eligible": False,
        "status": "passed",
    },
)
PY
module purge

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'GSE281364 Sei/MPRALegNet development evaluation: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
