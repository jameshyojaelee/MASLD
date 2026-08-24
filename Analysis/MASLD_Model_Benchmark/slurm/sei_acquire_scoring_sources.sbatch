#!/usr/bin/env bash
#SBATCH --job-name=masld-sei-scoring-src
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
SOURCE="${ROOT}/executions/alphagenome-sei-admission-sources-21064512"
SOURCE_SHA=84ae3470dca123c20904e315ef9926d5a139fb9fe6d99c51f6aaf694ef6e6785
REV=1e16bff9b1618af7065a23f82c4a3bb82aaeb689
TARGET="${ROOT}/executions/sei-scoring-sources-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
    printf 'Preserved failed Sei scoring-source stage: %s\n' "${FAILED}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite Sei scoring-source artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/sources" "${STAGE}/validation"
printf '%s  %s\n' "${SOURCE_SHA}" "${SOURCE}/ARTIFACTS.json" \
  | sha256sum --check --strict

BASE_URL="https://raw.githubusercontent.com/FunctionLab/sei-framework/${REV}"
for path in utils.py 2_raw_sc_score.py 2_varianteffect_sc_score.py; do
  curl --fail --location --retry 6 --retry-all-errors \
    --output "${STAGE}/sources/${path}" "${BASE_URL}/${path}"
done
printf '%s  %s\n' \
  d9d51305a65790c05643b82e4f41b13b6e9e67076893e30ef4b09a2d415a2b44 "${STAGE}/sources/utils.py" \
  0d0f7eefa02598a8ef4800ef36db0d65115f593a65a4c033fb3201698e71d141 "${STAGE}/sources/2_raw_sc_score.py" \
  df18dd5980bf56076798ee4acf95663620e867bf7889a896814a1f3b8035a630 "${STAGE}/sources/2_varianteffect_sc_score.py" \
  | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
python -m pip freeze > "${STAGE}/validation/effective-pip-freeze.txt"
PYTHONPATH="${ROOT}/src" python - "${SOURCE}" "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

source = Path(sys.argv[1]).resolve(strict=True)
root = Path(sys.argv[2]).resolve(strict=True)
verify_frozen_tree(source)
review = json.loads((source / "review/source_admission_receipt.json").read_text())
utils = (root / "sources/utils.py").read_text()
raw = (root / "sources/2_raw_sc_score.py").read_text()
variant = (root / "sources/2_varianteffect_sc_score.py").read_text()
if (
    review["status"] != "pass"
    or review["model_dispositions"]["sei"]["open_champion_eligible"]
    or "np.dot(chromatin_profile_preds, clustervfeat.T)" not in utils
    or "np.linalg.norm(clustervfeat, axis=1)" not in utils
    or "diffproj = altproj[:,:40] - refproj[:,:40]" not in utils
    or "sc_projection(input_preds, clustervfeat)" not in raw
    or "sc_hnorm_varianteffect(" not in variant
):
    raise SystemExit("Sei official scoring-source contract differs")
receipt = {
    "schema_version": "masld-bench-sei-scoring-source-v1",
    "status": "pass",
    "repository": "FunctionLab/sei-framework",
    "revision": "1e16bff9b1618af7065a23f82c4a3bb82aaeb689",
    "projection_array_shape": [61, 21907],
    "sequence_class_name_count": 40,
    "raw_projection_output_count": 61,
    "native_variant_output_count": 40,
    "native_variant_mapping": "first_40_projection_rows_after_histone_normalization_alt_minus_ref",
    "unlabeled_projection_rows": 21,
    "unlabeled_projection_rows_allowed_in_reported_native_variant_score": False,
    "terms": "academic_and_research_use_only_restricted_comparator",
    "model_or_checkpoint_loaded": False,
    "observed_outcomes_loaded": False,
    "sealed_outcomes_loaded": False,
}
(root / "validation/scoring_contract.json").write_text(
    json.dumps(receipt, sort_keys=True, indent=2) + "\n"
)
freeze_tree(
    root,
    {
        "artifact_class": "sei_official_scoring_sources",
        "upstream_source_artifacts_sha256": "84ae3470dca123c20904e315ef9926d5a139fb9fe6d99c51f6aaf694ef6e6785",
        "repository_revision": receipt["revision"],
        "native_variant_output_count": 40,
        "native_variant_mapping": receipt["native_variant_mapping"],
        "terms": receipt["terms"],
        "model_or_checkpoint_loaded": False,
    },
)
PY
module purge

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Sei official scoring sources: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
