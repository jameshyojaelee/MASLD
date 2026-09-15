#!/usr/bin/env bash
#SBATCH --job-name=model-data-054
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
export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=4

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
UNLABELED=${ROOT}/executions/model-data-051-21077382
UNLABELED_SHA256=ba40ce9555bb6d7bc278be2a80015b993f9f3650fa06f59505ad279c88cc117a
UNLABELED_H5_SHA256=b45f505fdbe737e7ad1332f47fd4394917c7b00002c95ad5b7172650caff59b5
SYMBOL_SOURCE=${ROOT}/executions/model-data-050-20260824T211318Z-21077049/fixture
SYMBOL_SOURCE_SHA256=672b1f75a1c3ca264066dd2a0bf24e479397d8ec5e2f67f13ccdd4be14879509
SYMBOL_H5_SHA256=6a704d4f0503d9bd113c9d70a6f92f65abd6af734a4f74cb2523c75c56d6a862
SPLIT=${ROOT}/executions/model-data-052-21077595
SPLIT_SHA256=f4a0fb205a8196aa01c07edfeb6d1db6833cf482be4b88e6f981a351e2fcd0e5
ROW_SPLIT_SHA256=f668e512ad116d945d58cd456e77e4e97472241a40c03f2e2c49b79e4477a819
MODEL=${ROOT}/executions/scimilarity-model-20260823/model_v1.1
GENE_ORDER_SHA256=0a9aeacf7fc19399426273d8cfab8683253ec56b51be556d42edd115f125b8ca
LIBRARY_MANIFEST=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/continual_integration/candidates/preflight-2026-08-13T224413Z-20206577/contract/library_manifest.tsv
LIBRARY_MANIFEST_SHA256=09dfc5668b512ac2f02fb876066d2e7c68be33f5c868de1521595f90809a8906
LIBRARY_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/integration/input_h5ad/human
TARGET=${ROOT}/executions/scimilarity-v1-1-frozen-screen-fixture-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed
PYTHON=${RUNTIME}/env/bin/python

finish() {
  exit_code=$?
  if [[ "${exit_code}" -ne 0 && -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
  exit "${exit_code}"
}
trap finish EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite SCimilarity frozen-screen fixture.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"
printf '%s  %s\n' \
  "${RUNTIME_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  "${UNLABELED_SHA256}" "${UNLABELED}/ARTIFACTS.json" \
  "${UNLABELED_H5_SHA256}" "${UNLABELED}/fixture/unlabeled_counts.h5ad" \
  "${SYMBOL_SOURCE_SHA256}" "${SYMBOL_SOURCE}/ARTIFACTS.json" \
  "${SYMBOL_H5_SHA256}" "${SYMBOL_SOURCE}/resource_atlas_frozen_screen_50000.h5ad" \
  "${SPLIT_SHA256}" "${SPLIT}/ARTIFACTS.json" \
  "${ROW_SPLIT_SHA256}" "${SPLIT}/split/row_outer_folds.tsv" \
  "${GENE_ORDER_SHA256}" "${MODEL}/gene_order.tsv" \
  "${LIBRARY_MANIFEST_SHA256}" "${LIBRARY_MANIFEST}" \
  > "${STAGE}/validation/input.sha256"
sha256sum --check --strict "${STAGE}/validation/input.sha256"
export PYTHONPATH="${ROOT}:${ROOT}/src"
"${PYTHON}" -m unittest -v tests.unit.test_scimilarity_v1_1_frozen_screen_50000 \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
"${PYTHON}" scripts/build_scimilarity_v1_1_frozen_screen_50000.py \
  --source "${UNLABELED}/fixture/unlabeled_counts.h5ad" \
  --symbol-source "${SYMBOL_SOURCE}/resource_atlas_frozen_screen_50000.h5ad" \
  --split "${SPLIT}/split/row_outer_folds.tsv" \
  --gene-order "${MODEL}/gene_order.tsv" \
  --library-manifest "${LIBRARY_MANIFEST}" \
  --library-root "${LIBRARY_ROOT}" \
  --output "${STAGE}/fixture" \
  --expected-source-sha256 "${UNLABELED_H5_SHA256}" \
  --expected-symbol-source-sha256 "${SYMBOL_H5_SHA256}" \
  --expected-split-sha256 "${ROW_SPLIT_SHA256}" \
  --expected-gene-order-sha256 "${GENE_ORDER_SHA256}" \
  --expected-library-manifest-sha256 "${LIBRARY_MANIFEST_SHA256}" \
  > "${STAGE}/validation/fixture.stdout.json" \
  2> "${STAGE}/validation/fixture.stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/validation/pip-freeze.txt"
sha256sum \
  scripts/build_scimilarity_v1_1_frozen_screen_50000.py \
  scripts/extract_scimilarity_v1_1_frozen_screen_50000.py \
  tests/unit/test_scimilarity_v1_1_frozen_screen_50000.py \
  slurm/build_scimilarity_v1_1_frozen_screen_50000_cpu.sbatch \
  "${RUNTIME}/ARTIFACTS.json" \
  "${UNLABELED}/ARTIFACTS.json" \
  "${SYMBOL_SOURCE}/ARTIFACTS.json" \
  "${SPLIT}/ARTIFACTS.json" \
  > "${STAGE}/validation/source.sha256"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "fixture/fixture_receipt.json").read_text(encoding="utf-8"))
if (
    receipt.get("status") != "pass_outcome_blind_fixture"
    or receipt.get("rows") != 50_000
    or receipt.get("donors") != 102
    or receipt.get("studies") != 7
    or receipt.get("input_dimension") != 28_231
    or receipt.get("gene_overlap", 0) < 5_000
    or receipt.get("evaluation_label_columns_read")
    or receipt.get("histology_columns_read")
    or receipt.get("sealed_outcomes_read") is not False
    or receipt.get("released_reference_index_used") is not False
    or receipt.get("aggregate_exposure_status") != "encoder_seen"
    or receipt.get("study_gene_observability_shape") != [7, 28_231]
    or receipt.get("native_model_accepts_observability_mask") is not False
    or receipt.get("feature_symbol_source_observations_read") is not False
    or receipt.get("exposure_by_study", {}).get("GSE185477") != "encoder_seen"
    or receipt.get("exposure_by_study", {}).get("GSE136103") != "reference_only"
):
    raise SystemExit("SCimilarity 50k fixture receipt differs")
freeze_tree(
    root,
    {
        "artifact_class": "scimilarity_v1_1_frozen_screen_input",
        "dataset_view_id": "resource_atlas_frozen_screen_50000_v1",
        "development_rows": 50_000,
        "donors": 102,
        "studies": 7,
        "evaluation_labels_read": False,
        "feature_symbol_source_observations_read": False,
        "released_reference_index_used": False,
        "aggregate_exposure_status": "encoder_seen",
        "sealed_outcomes_read": False,
        "status": "passed",
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
