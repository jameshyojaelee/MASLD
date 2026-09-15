#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-train-003
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4

: "${PSEUDOBULK:?PSEUDOBULK path is required}"
: "${PSEUDOBULK_SHA256:?PSEUDOBULK_SHA256 is required}"

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PYTHON="${PROJECT}/.agent/scanpy_env/bin/python"
SPEC="${ROOT}/executions/gse296875-phenotype-campaign-spec-21107179"
TARGET="${ROOT}/executions/gse296875-phenotype-predictions-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi
}
trap preserve_failed_stage EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite prediction artifact: %s\n' "${TARGET}" >&2
  exit 1
fi

cd "${ROOT}"
mkdir -p "${STAGE}/validation"

cat > "${STAGE}/validation/bound_fixtures.sha256" <<SUMS
6518494a7249cb1eac81fbad75989e7dfe0269b611254d5097736b9f5b1da140  executions/gse296875-phenotype-endpoints-20260825/ARTIFACTS.json
e6bc9161acaab508c298c49036af1ad10a720949178a14cba42282933195d416  executions/gse296875-donor-folds-20260825/ARTIFACTS.json
474db3b44b22ff3c3556d1fcca92c5c17d8a10f2b873b5c7ae345c81d6382a58  executions/gse296875-phenotype-campaign-spec-21107179/ARTIFACTS.json
${PSEUDOBULK_SHA256}  ${PSEUDOBULK}/ARTIFACTS.json
SUMS
sha256sum --check --strict "${STAGE}/validation/bound_fixtures.sha256"

"${PYTHON}" -m unittest -v tests.unit.test_fit_gse296875_phenotype_predictions \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

"${PYTHON}" scripts/fit_gse296875_phenotype_predictions.py \
  --root "${ROOT}" \
  --pseudobulk "${PSEUDOBULK}" \
  --pseudobulk-sha256 "${PSEUDOBULK_SHA256}" \
  --output "${STAGE}/predictions" \
  > "${STAGE}/validation/fit_stdout.json"

"${PYTHON}" -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt" 2>&1 || true

sha256sum \
  scripts/fit_gse296875_phenotype_predictions.py \
  tests/unit/test_fit_gse296875_phenotype_predictions.py \
  slurm/fit_gse296875_phenotype_predictions_cpu.sbatch \
  > "${STAGE}/source.sha256"

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH=src python3 - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
audit = json.loads((root / "predictions" / "audit.json").read_text())
if audit["metrics_computed"] != 0:
    raise SystemExit("the fitting process must not compute a metric")
if audit["test_fold_labels_never_entered_a_fit_or_a_transform"] is not True:
    raise SystemExit("fold firewall assertion missing")
if audit["masked_units_never_zero"] is not True:
    raise SystemExit("masked units must never be encoded as zero")
freeze_tree(root, {
    "artifact_class": "donor_lineage_phenotype_out_of_fold_predictions",
    "dataset_id": "gse296875",
    "donors": 39,
    "unit_of_inference": "donor",
    "campaign_spec_sha256": "474db3b44b22ff3c3556d1fcca92c5c17d8a10f2b873b5c7ae345c81d6382a58",
    "metrics_computed": 0,
    "frozen_before_any_outcome_join": True,
    "status": "passed",
})
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'GSE296875 phenotype predictions: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
