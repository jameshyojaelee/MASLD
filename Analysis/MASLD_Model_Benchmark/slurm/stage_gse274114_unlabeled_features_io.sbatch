#!/usr/bin/env bash
#SBATCH --job-name=data-stage-242
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --account=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=executions/logs/data-stage-242-%j.out
#SBATCH --error=executions/logs/data-stage-242-%j.err

set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
VALIDATION=${ROOT}/executions/data-check-241-21092012
VALIDATION_SHA256=40d2d58f94344e86a31764a229d66090836f58046f5536a090637e1645eb3737
FIXTURE=${ROOT}/executions/model-data-080-21082249-gse274114
FIXTURE_SHA256=b334ff4fe77903be82bb228fc479a28d92d74f4cec0b06ea98c3d6a6100c9300
TARGET=${ROOT}/executions/data-stage-242-${SLURM_JOB_ID}-gse274114-unlabeled
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d ${STAGE} && ! -e ${FAILED} ]]; then mv "${STAGE}" "${FAILED}"; fi
}
trap preserve_failed_stage EXIT
if [[ -e ${TARGET} || -e ${STAGE} || -e ${FAILED} ]]; then
  printf 'Refusing to overwrite GSE274114 unlabeled feature staging.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
cd "${ROOT}"
umask 027
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 PYTHONPATH="${ROOT}/src:${ROOT}"
printf '%s  %s\n' "${RUNTIME_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${VALIDATION_SHA256}" "${VALIDATION}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${FIXTURE_SHA256}" "${FIXTURE}/ARTIFACTS.json" \
  | sha256sum --check --strict
sha256sum --check --strict "${VALIDATION}/source.sha256"

"${PYTHON}" scripts/gse274114_baseline_firewall.py stage-features \
  --fixture-root "${FIXTURE}" \
  --output "${STAGE}/features" \
  > "${STAGE}/stage.stdout.json" 2> "${STAGE}/stage.stderr.txt"
"${PYTHON}" - "${STAGE}/features" <<'PY'
from pathlib import Path
import csv
import json
import sys

import numpy as np

from masld_bench.artifacts import verify_frozen_tree

root = Path(sys.argv[1])
manifest = verify_frozen_tree(root)
receipt = json.loads((root / "receipt.json").read_text(encoding="utf-8"))
expected = {
    "participants": 39,
    "genes": 60324,
    "labels_opened": False,
    "folds_opened": False,
    "participant_metadata_opened": False,
    "molecular_values_opened": True,
    "normalization_run": False,
    "model_fit_or_scoring_run": False,
    "effective_length_retained": True,
}
if (
    manifest["metadata"].get("artifact_class")
    != "gse274114_unlabeled_feature_matrix"
    or any(receipt.get(key) != value for key, value in expected.items())
    or any(manifest["metadata"].get(key) != value for key, value in expected.items())
):
    raise SystemExit("unlabeled feature-stage firewall differs")
with (root / "participants.tsv").open(encoding="utf-8", newline="") as handle:
    participant_rows = list(csv.DictReader(handle, delimiter="\t"))
with (root / "genes.tsv").open(encoding="utf-8", newline="") as handle:
    gene_rows = list(csv.DictReader(handle, delimiter="\t"))
counts = np.load(root / "numreads.npy", mmap_mode="r", allow_pickle=False)
lengths = np.load(root / "effective_length.npy", mmap_mode="r", allow_pickle=False)
if (
    len(participant_rows) != 39
    or tuple(participant_rows[0]) != ("row_id",)
    or len(gene_rows) != 60324
    or tuple(gene_rows[0]) != ("gene_id",)
    or counts.shape != (39, 60324)
    or lengths.shape != counts.shape
    or counts.dtype != np.float64
    or lengths.dtype != np.float64
):
    raise SystemExit("unlabeled feature-stage dimensions differ")
PY
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum \
  scripts/gse274114_baseline_firewall.py \
  config/campaigns/gse274114_within_platform_baselines_v1.toml \
  slurm/stage_gse274114_unlabeled_features_io.sbatch \
  "${VALIDATION}/ARTIFACTS.json" \
  "${FIXTURE}/ARTIFACTS.json" \
  "${STAGE}/features/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1])
verify_frozen_tree(root / "features")
receipt = json.loads((root / "features/receipt.json").read_text(encoding="utf-8"))
freeze_tree(
    root,
    {
        "artifact_class": "gse274114_outcome_blind_unlabeled_feature_stage",
        "status": "passed_outcome_blind_molecular_stage",
        "participants": receipt["participants"],
        "genes": receipt["genes"],
        "real_evaluator_labels_opened": False,
        "real_evaluator_folds_opened": False,
        "real_participant_metadata_opened": False,
        "real_molecular_values_opened": True,
        "normalization_run": False,
        "model_fit_or_scoring_run": False,
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/features/ARTIFACTS.json"
