#!/usr/bin/env bash
#SBATCH --job-name=model-data-051
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
SOURCE_ROOT=${ROOT}/executions/model-data-050-20260824T211318Z-21077049/fixture
SOURCE_SHA256=6a704d4f0503d9bd113c9d70a6f92f65abd6af734a4f74cb2523c75c56d6a862
SOURCE_ARTIFACTS_SHA256=672b1f75a1c3ca264066dd2a0bf24e479397d8ec5e2f67f13ccdd4be14879509
SOURCE=${SOURCE_ROOT}/resource_atlas_frozen_screen_50000.h5ad
CHECKPOINT=${ROOT}/executions/transcriptformer-materialized-21066178/bundle/tf_sapiens
VOCAB=${CHECKPOINT}/vocabs/homo_sapiens_gene.h5
TARGET=${ROOT}/executions/model-data-051-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed
PYTHON=${RUNTIME}/env/bin/python

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite frozen-screen token artifacts.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
printf '%s  %s\n' "${RUNTIME_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${SOURCE_ARTIFACTS_SHA256}" "${SOURCE_ROOT}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${SOURCE_SHA256}" "${SOURCE}" \
  | sha256sum --check --strict

cd "${ROOT}"
export PYTHONPATH="${ROOT}:${ROOT}/src"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
"${PYTHON}" -m unittest tests.unit.test_build_transcriptformer_unlabeled_fixture \
  > "${STAGE}/validation/unit_tests.txt" 2>&1
"${PYTHON}" - "${SOURCE_ROOT}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree

verified = verify_frozen_tree(Path(sys.argv[1]))
if verified["metadata"].get("row_count") != 50_000:
    raise SystemExit("frozen-screen source row count differs")
PY
"${PYTHON}" scripts/build_transcriptformer_unlabeled_fixture.py \
  --source "${SOURCE}" \
  --vocab "${VOCAB}" \
  --output "${STAGE}/fixture" \
  --expected-source-sha256 "${SOURCE_SHA256}" \
  --expected-rows 50000 \
  --expected-donors 102 \
  --dataset-view-id resource_atlas_frozen_screen_50000_v1 \
  --fold-namespace resource_atlas_current:cell_state_mapping:donor_outer:v1 \
  > "${STAGE}/validation/fixture_stdout.json" \
  2> "${STAGE}/validation/fixture_stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt"
sha256sum \
  scripts/build_transcriptformer_unlabeled_fixture.py \
  tests/unit/test_build_transcriptformer_unlabeled_fixture.py \
  slurm/build_transcriptformer_frozen_screen_50000.sbatch \
  "${RUNTIME}/ARTIFACTS.json" \
  "${SOURCE_ROOT}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

import h5py
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
fixture = root / "fixture"
receipt = json.loads((fixture / "fixture_receipt.json").read_text(encoding="utf-8"))
if (
    receipt.get("status") != "pass_outcome_blind_fixture"
    or receipt.get("rows") != 50_000
    or receipt.get("donors") != 102
    or receipt.get("dataset_view_id") != "resource_atlas_frozen_screen_50000_v1"
    or receipt.get("outer_fold_namespace")
    != "resource_atlas_current:cell_state_mapping:donor_outer:v1"
    or receipt.get("evaluation_label_columns_read")
    or receipt.get("sealed_outcomes_read")
    or receipt.get("vocabulary_rows_in_checkpoint") != 23_830
):
    raise SystemExit("TranscriptFormer 50k fixture receipt differs")
with h5py.File(fixture / "unlabeled_counts.h5ad", "r") as handle:
    fields = set(handle["obs"].keys())
    index_key = handle["obs"].attrs.get("_index", "_index")
    if isinstance(index_key, bytes):
        index_key = index_key.decode("utf-8")
allowed = {str(index_key), "row_id", "donor_id", "dataset_id", "outer_fold", "assay"}
if not fields.issubset(allowed) or not {"donor_id", "dataset_id", "outer_fold"}.issubset(fields):
    raise SystemExit("TranscriptFormer 50k label firewall differs")
freeze_tree(root, {
    "artifact_class": "transcriptformer_frozen_screen_input",
    "dataset_view_id": "resource_atlas_frozen_screen_50000_v1",
    "development_rows": 50_000,
    "donors": 102,
    "evaluation_labels_read": False,
    "fold_namespace": "resource_atlas_current:cell_state_mapping:donor_outer:v1",
    "gene_order_policies": ["native", "common"],
    "sealed_outcomes_read": False,
    "status": "passed",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
