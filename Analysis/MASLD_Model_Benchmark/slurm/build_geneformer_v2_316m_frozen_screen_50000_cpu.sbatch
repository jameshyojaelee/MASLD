#!/usr/bin/env bash
#SBATCH --job-name=model-data-053
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
export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1
export WANDB_MODE=offline
export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=4

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/geneformer-v2-316m-compat-runtime-21077474
RUNTIME_SHA256=c034ee71fb2c67f04f11e91b41ea8c9bf6a7b9df3a4bb72378e5b2696058c307
RUNTIME_ARCHIVE_SHA256=a4fdc1f1d0b67c4753354c5a04ccc05c84a3bf7442d8655bc8ed6493a2487775
RESOURCE=${ROOT}/executions/model-data-050-20260824T211318Z-21077049/fixture
RESOURCE_SHA256=672b1f75a1c3ca264066dd2a0bf24e479397d8ec5e2f67f13ccdd4be14879509
RESOURCE_H5_SHA256=6a704d4f0503d9bd113c9d70a6f92f65abd6af734a4f74cb2523c75c56d6a862
UNLABELED=${ROOT}/executions/model-data-051-21077382
UNLABELED_SHA256=ba40ce9555bb6d7bc278be2a80015b993f9f3650fa06f59505ad279c88cc117a
UNLABELED_H5_SHA256=b45f505fdbe737e7ad1332f47fd4394917c7b00002c95ad5b7172650caff59b5
SPLIT=${ROOT}/executions/model-data-052-21077595
SPLIT_SHA256=f4a0fb205a8196aa01c07edfeb6d1db6833cf482be4b88e6f981a351e2fcd0e5
ROW_SPLIT_SHA256=f668e512ad116d945d58cd456e77e4e97472241a40c03f2e2c49b79e4477a819
BUNDLE=${ROOT}/executions/geneformer-acquisition-21014827/bundles/geneformer_v2_316m--11f311a5078fc4c8
BUNDLE_SHA256=b269cef53fa92cc7e7ab9c21a2b85de69904cdfaaf5881e9676c57dc2ae02ae0
TARGET=${ROOT}/executions/geneformer-v2-316m-frozen-screen-fixture-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed
WORK=$(mktemp -d "${TMPDIR:-/tmp}/geneformer-50k-fixture.${SLURM_JOB_ID}.XXXXXX")

finish() {
  exit_code=$?
  if [[ "${exit_code}" -ne 0 && -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
  if [[ -d "${WORK}" ]]; then
    rm -rf -- "${WORK}"
  fi
  exit "${exit_code}"
}
trap finish EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite Geneformer 50k fixture.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation" "${WORK}/runtime"
cd "${ROOT}"

printf '%s  %s\n' \
  "${RUNTIME_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  "${RUNTIME_ARCHIVE_SHA256}" "${RUNTIME}/payload/geneformer-decima-compat-env.tar.gz" \
  "${RESOURCE_SHA256}" "${RESOURCE}/ARTIFACTS.json" \
  "${RESOURCE_H5_SHA256}" "${RESOURCE}/resource_atlas_frozen_screen_50000.h5ad" \
  "${UNLABELED_SHA256}" "${UNLABELED}/ARTIFACTS.json" \
  "${UNLABELED_H5_SHA256}" "${UNLABELED}/fixture/unlabeled_counts.h5ad" \
  "${SPLIT_SHA256}" "${SPLIT}/ARTIFACTS.json" \
  "${ROW_SPLIT_SHA256}" "${SPLIT}/split/row_outer_folds.tsv" \
  "${BUNDLE_SHA256}" "${BUNDLE}/ARTIFACTS.json" \
  > "${STAGE}/validation/input.sha256"
sha256sum --check --strict "${STAGE}/validation/input.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${RESOURCE}" "${UNLABELED}" "${SPLIT}" "${BUNDLE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

resource, unlabeled, split, bundle = map(Path, sys.argv[1:])
for path in (resource, unlabeled, split, bundle):
    verify_frozen_tree(path)
PY
module purge

tar -xzf "${RUNTIME}/payload/geneformer-decima-compat-env.tar.gz" -C "${WORK}/runtime"
export PYTHONPATH="${ROOT}:${ROOT}/src"
"${WORK}/runtime/bin/python" -m unittest -v \
  tests.unit.test_geneformer_v2_316m_frozen_screen_50000 \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
"${WORK}/runtime/bin/python" \
  scripts/build_geneformer_v2_316m_frozen_screen_50000.py \
  --source "${UNLABELED}/fixture/unlabeled_counts.h5ad" \
  --split "${SPLIT}/split/row_outer_folds.tsv" \
  --median "${BUNDLE}/sanitized/geneformer/gene_median_dictionary_gc104M.pkl.json" \
  --token-dictionary "${BUNDLE}/sanitized/geneformer/token_dictionary_gc104M.pkl.json" \
  --config "${BUNDLE}/upstream/Geneformer-V2-316M/config.json" \
  --output "${STAGE}/fixture" \
  --expected-source-sha256 "${UNLABELED_H5_SHA256}" \
  --expected-split-sha256 "${ROW_SPLIT_SHA256}" \
  > "${STAGE}/validation/builder.stdout.txt" \
  2> "${STAGE}/validation/builder.stderr.txt"
"${WORK}/runtime/bin/python" -m pip freeze --all \
  > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  scripts/build_geneformer_v2_316m_frozen_screen_50000.py \
  scripts/extract_geneformer_v2_316m_frozen_screen_50000.py \
  tests/unit/test_geneformer_v2_316m_frozen_screen_50000.py \
  slurm/build_geneformer_v2_316m_frozen_screen_50000_cpu.sbatch \
  "${RUNTIME}/ARTIFACTS.json" \
  "${RESOURCE}/ARTIFACTS.json" \
  "${UNLABELED}/ARTIFACTS.json" \
  "${SPLIT}/ARTIFACTS.json" \
  "${BUNDLE}/ARTIFACTS.json" \
  > "${STAGE}/validation/source.sha256"

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

stage = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((stage / "fixture/fixture_receipt.json").read_text(encoding="utf-8"))
if (
    receipt.get("status") != "pass_outcome_blind_fixture"
    or receipt.get("rows") != 50_000
    or receipt.get("donors") != 102
    or receipt.get("studies") != 7
    or receipt.get("split_id") != "resource_atlas_study_outer_5fold_v1"
    or receipt.get("evaluation_label_columns_read")
    or receipt.get("histology_columns_read")
    or receipt.get("sealed_outcomes_read") is not False
):
    raise SystemExit("Geneformer 50k fixture receipt differs")
freeze_tree(
    stage,
    {
        "artifact_class": "geneformer_v2_316m_frozen_screen_fixture",
        "model_id": "geneformer_v2_316m",
        "dataset_view_id": "resource_atlas_frozen_screen_50000_v1",
        "split_id": "resource_atlas_study_outer_5fold_v1",
        "rows": 50_000,
        "donors": 102,
        "studies": 7,
        "outcome_blind": True,
        "evaluation_labels_read": False,
        "histology_read": False,
        "sealed_outcomes_read": False,
        "status": "passed",
    },
)
PY

trap - EXIT
rm -rf -- "${WORK}"
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
