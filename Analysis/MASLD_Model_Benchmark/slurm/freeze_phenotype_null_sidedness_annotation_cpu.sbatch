#!/usr/bin/env bash
#SBATCH --job-name=model-check-013
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TARGET="${ROOT}/executions/phenotype-null-sidedness-annotation-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"; FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite annotation artifact: %s\n' "${TARGET}" >&2; exit 1
fi
cd "${ROOT}"
module load python/3.11.5-GCCcore-13.2.0
mkdir -p "${STAGE}/validation" "${STAGE}/declarations"

PYTHONPATH=src python3 -m unittest -v \
  tests.unit.test_ranking_metric_floor_requirement \
  tests.unit.test_gse296875_phenotype_scoring \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

cp config/campaigns/phenotype_null_sidedness_annotation_20260825.json "${STAGE}/declarations/"
python3 -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt" 2>&1 || true
sha256sum \
  config/campaigns/phenotype_null_sidedness_annotation_20260825.json \
  src/masld_bench/gse296875_phenotype_scoring.py \
  src/masld_bench/evaluators/auprc_reference.py \
  tests/unit/test_ranking_metric_floor_requirement.py \
  slurm/freeze_phenotype_null_sidedness_annotation_cpu.sbatch \
  > "${STAGE}/source.sha256"

PYTHONPATH=src python3 - "${STAGE}" <<'PY'
from pathlib import Path
import json, sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
spec = json.loads((root / "declarations" / "phenotype_null_sidedness_annotation_20260825.json").read_text())
if spec["record_only"] is not True or spec["consumed_by_runtime"] is not False:
    raise SystemExit("the annotation must be record-only")
if spec["reported_metrics_unchanged"] is not True:
    raise SystemExit("the annotation must change no reported metric")
if spec["impact_on_conclusions"]["any_conclusion_changes"] is not False:
    raise SystemExit("annotation claims a conclusion changed")
freeze_tree(root, {
    "artifact_class": "null_sidedness_record_only_annotation",
    "dataset_id": "gse296875",
    "record_only": True,
    "reported_metrics_unchanged": True,
    "status": "passed",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Null sidedness annotation: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
