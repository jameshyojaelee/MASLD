#!/usr/bin/env bash
#SBATCH --job-name=model-data-064
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
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
ACTIVATION_CAMPAIGN=${ROOT}/executions/model-data-061-21079623
ACTIVATION_VIEW=${ACTIVATION_CAMPAIGN}/activation
MATRIX=${ROOT}/executions/gse267145-matrix-audit-21064333
JOIN=${ROOT}/executions/gse267145-authoritative-join-21064930
MEASUREMENT=${ROOT}/executions/gse267145-rna-measurement-21065097
REFERENCE=${ROOT}/executions/gse267145-reference-crosswalk-21065154
QC_RIGHTS=${ROOT}/executions/gse267145-qc-rights-21065336
COORDINATE=${ROOT}/executions/gse267145-coordinate-activation-21070593
TARGET=${ROOT}/executions/model-data-064-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite paired molecular fixture artifacts.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
cd "${ROOT}"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}/src:${ROOT}"
printf '%s  %s\n' "${RUNTIME_ARTIFACTS_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict

"${PYTHON}" -m unittest -v tests.unit.test_build_gse267145_paired_molecular_fixture \
  > "${STAGE}/unit_tests.txt" 2>&1
"${PYTHON}" - "${ACTIVATION_CAMPAIGN}" "${ACTIVATION_VIEW}" "${MATRIX}" "${JOIN}" \
  "${MEASUREMENT}" "${REFERENCE}" "${QC_RIGHTS}" "${COORDINATE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

for value in sys.argv[1:]:
    verify_frozen_tree(Path(value).resolve(strict=True))
print("all upstream frozen-tree verifications passed")
PY
"${PYTHON}" scripts/build_gse267145_paired_molecular_fixture.py \
  --activation-campaign "${ACTIVATION_CAMPAIGN}" \
  --activation-campaign-artifacts-sha256 e62c882baaaa236af2a6d4a6818abc458ad0810c30014969ccee9d055b8df299 \
  --activation-view "${ACTIVATION_VIEW}" \
  --activation-view-artifacts-sha256 98a74b97f6a9712a109d2516c3803a4acf15030ed607f861e1ce8c22b4f9753a \
  --matrix-audit "${MATRIX}" \
  --matrix-artifacts-sha256 4e0d17f392530f4e6f0b54346d5d9333e450b7606d9698fbf7253a4f2a0e0417 \
  --join "${JOIN}" \
  --join-artifacts-sha256 e6c539c5fb29cd357dc126073779948c73c27fd219ff7a4d80b3d20a6cddb580 \
  --measurement "${MEASUREMENT}" \
  --measurement-artifacts-sha256 38531e463837cafd2c8e321411048c7b8b493876c3c964c0f50790cac8806f2f \
  --reference "${REFERENCE}" \
  --reference-artifacts-sha256 da491302ea738a92017978518239bc3884a7635d580e01e3bd28cd4f55b630c2 \
  --qc-rights "${QC_RIGHTS}" \
  --qc-rights-artifacts-sha256 ac7c92fd01b2b4240b380e4f6a58d98f1a6f853f8d0ad98d7167d7af34baaa73 \
  --coordinate-activation "${COORDINATE}" \
  --coordinate-artifacts-sha256 e5042eb43dddd84cda81c5aeed5d6dbeb0378ab1c942a0183e2ef1b41cd18b20 \
  --output "${STAGE}/fixture" \
  > "${STAGE}/run_stdout.json" 2> "${STAGE}/run_stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum \
  scripts/build_gse267145_paired_molecular_fixture.py \
  tests/unit/test_build_gse267145_paired_molecular_fixture.py \
  slurm/build_gse267145_paired_molecular_fixture_cpu.sbatch \
  "${RUNTIME}/ARTIFACTS.json" \
  "${ACTIVATION_CAMPAIGN}/ARTIFACTS.json" \
  "${ACTIVATION_VIEW}/ARTIFACTS.json" \
  "${MATRIX}/ARTIFACTS.json" \
  "${JOIN}/ARTIFACTS.json" \
  "${MEASUREMENT}/ARTIFACTS.json" \
  "${REFERENCE}/ARTIFACTS.json" \
  "${QC_RIGHTS}/ARTIFACTS.json" \
  "${COORDINATE}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import csv
import json
import sys

import numpy as np

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
fixture = root / "fixture"
molecular = fixture / "molecular"
folds = fixture / "folds"
receipt = json.loads((molecular / "receipt.json").read_text(encoding="utf-8"))
rna = np.load(molecular / "rna_values.npy", mmap_mode="r", allow_pickle=False)
h3 = np.load(molecular / "h3k27ac_counts.npy", mmap_mode="r", allow_pickle=False)
rna_mask = np.load(molecular / "rna_observed_mask.npy", allow_pickle=False)
h3_mask = np.load(molecular / "h3k27ac_observed_mask.npy", allow_pickle=False)
with (molecular / "participant_axis.tsv").open(encoding="utf-8", newline="") as handle:
    participant_rows = list(csv.DictReader(handle, delimiter="\t"))
with (folds / "participant_outer_folds.tsv").open(encoding="utf-8", newline="") as handle:
    fold_rows = list(csv.DictReader(handle, delimiter="\t"))
if (
    rna.shape != (99, 42163)
    or rna.dtype != np.float64
    or not np.any(rna != np.floor(rna))
    or h3.shape != (99, 96460)
    or h3.dtype != np.uint32
    or rna_mask.dtype != np.bool_
    or h3_mask.dtype != np.bool_
    or not rna_mask.all()
    or not h3_mask.all()
    or len(participant_rows) != 99
    or len(fold_rows) != 99
    or [row["participant_id"] for row in participant_rows]
    != [row["participant_id"] for row in fold_rows]
    or receipt.get("feature_firewall", {}).get("outcome_values_accessed") is not False
    or receipt.get("model_fitted") is not False
    or receipt.get("metrics_calculated") is not False
):
    raise SystemExit("paired molecular fixture verification differs")
freeze_tree(
    molecular,
    {
        "artifact_class": "gse267145_outcome_free_paired_molecular_values",
        "participants": 99,
        "rna_features": 42163,
        "h3k27ac_opaque_features": 96460,
        "pairing": "same_sample_different_aliquot",
        "explicit_modality_masks": True,
        "outcome_values_included": False,
        "model_fitted": False,
        "metrics_calculated": False,
        "status": "passed",
    },
)
verify_frozen_tree(molecular)
freeze_tree(
    folds,
    {
        "artifact_class": "gse267145_participant_outer_folds",
        "participants": 99,
        "biological_unit": "participant",
        "folds_recomputed": False,
        "outcome_values_included": False,
        "status": "passed",
    },
)
verify_frozen_tree(folds)
freeze_tree(
    fixture,
    {
        "artifact_class": "gse267145_outcome_free_paired_molecular_fixture",
        "participants": 99,
        "rna_features": 42163,
        "h3k27ac_opaque_features": 96460,
        "molecular_and_folds_frozen_separately": True,
        "outcome_values_included": False,
        "model_fitted": False,
        "metrics_calculated": False,
        "status": "passed",
    },
)
verify_frozen_tree(fixture)
freeze_tree(
    root,
    {
        "artifact_class": "gse267145_paired_molecular_fixture_campaign",
        "participants": 99,
        "model_fitted": False,
        "metrics_calculated": False,
        "status": "passed",
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum \
  "${TARGET}/ARTIFACTS.json" \
  "${TARGET}/fixture/ARTIFACTS.json" \
  "${TARGET}/fixture/molecular/ARTIFACTS.json" \
  "${TARGET}/fixture/folds/ARTIFACTS.json"
