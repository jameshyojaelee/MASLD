#!/usr/bin/env bash
#SBATCH --job-name=model-work-202
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export CUBLAS_WORKSPACE_CONFIG=:4096:8
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
BASE="${ROOT}/executions/environments/scvi-tools-1.5.0.post1-portable-21064579"
TORCH="${ROOT}/executions/environments/pytorch-2.6.0-cu124-portable-21064747"
CONVERSION="${ROOT}/executions/environments/dna-lm-torch2.6-conversion-21064945"
SOURCE="${ROOT}/executions/alphagenome-sei-admission-sources-21064512"
SCORING="${ROOT}/executions/sei-scoring-sources-21065483"
WEIGHTS="${ROOT}/executions/sei-weights-only-conversion-21065468"
FIXTURE="${ROOT}/executions/gse281364-sei-fixture-21068925"
RUNTIME_PROBE="${ROOT}/executions/sei-native-probe-21066123"
OFFLINE_GUARD="${RUNTIME_PROBE}/offline-guard"
TARGET="${ROOT}/executions/sei-gse281364-smoke-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite Sei prediction artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
for pair in \
  "30814d4944e933d999e8c02c5b521432786b1985cc8794f9436f2a44af5ad1af ${BASE}/ARTIFACTS.json" \
  "d4112dbc9574e619edbaec0f2c1503c57e0271fe17517b2e121d711ea21f71a3 ${TORCH}/ARTIFACTS.json" \
  "5fed22c975200495ea53b38b09f5286cf13a9230a138cb72a0c758c85d5cb9d4 ${CONVERSION}/ARTIFACTS.json" \
  "84ae3470dca123c20904e315ef9926d5a139fb9fe6d99c51f6aaf694ef6e6785 ${SOURCE}/ARTIFACTS.json" \
  "529a2034a4a8168d8e436c7044364c65735bb274a43b9652d65000110562c5f7 ${SCORING}/ARTIFACTS.json" \
  "2b19a580f496db0a30a6174cc4f7fe65b0cf399e5dde94c8087665d374ba892b ${WEIGHTS}/ARTIFACTS.json" \
  "620464c0ef7cf378eb2b44a1608cb6b9c273bd70b8770a43b2a95d21e725af4d ${FIXTURE}/ARTIFACTS.json" \
  "c0607cba297ad2e3c2889ebd4825b287655e6f0064f0d9137aaa60be34026894 ${RUNTIME_PROBE}/ARTIFACTS.json"; do
  read -r expected path <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${path}" | sha256sum --check --strict
done
for pair in \
  "37b83ff34ef07fb724098bb229d532903af81008f760e934795dc5de77c31d70 ${ROOT}/scripts/sei_predict_gse281364.py" \
  "f124a901653446e67a5914d96f0e8fd9e3de013b578eef5887b25220267899ce ${ROOT}/tests/unit/test_sei_predict_gse281364.py"; do
  read -r expected path <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${path}" | sha256sum --check --strict
done
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${BASE}" "${TORCH}" "${CONVERSION}" "${SOURCE}" "${SCORING}" "${WEIGHTS}" "${FIXTURE}" "${RUNTIME_PROBE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for value in sys.argv[1:]:
    verify_frozen_tree(Path(value).resolve(strict=True))
PY
if [[ "${PREFLIGHT_ONLY:-0}" == 1 ]]; then
  printf 'Sei GSE281364 prediction preflight passed\n'
  exit 0
fi
install -d -m 0750 "${STAGE}/validation"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${OFFLINE_GUARD}:${CONVERSION}/site-packages:${TORCH}/site-packages:${BASE}/site-packages:${ROOT}"
python -m unittest -v tests.unit.test_sei_predict_gse281364 \
  > "${STAGE}/validation/tests.stdout.txt" 2> "${STAGE}/validation/tests.stderr.txt"
python "${ROOT}/scripts/sei_predict_gse281364.py" \
  --model-source "${SOURCE}/sources/sei_model_py" \
  --scoring-source "${SCORING}/sources/utils.py" \
  --checkpoint "${WEIGHTS}/converted/sei.model.safetensors" \
  --projection "${WEIGHTS}/native/projvec_targets.npy" \
  --histone-indices "${WEIGHTS}/native/histone_inds.npy" \
  --fixture-artifacts "${FIXTURE}/ARTIFACTS.json" \
  --manifest "${FIXTURE}/fixture/manifest.tsv" \
  --fasta "${FIXTURE}/fixture/sei.alleles.fa.gz" \
  --batch-variants 4 --output "${STAGE}/predictions" \
  > "${STAGE}/validation/prediction.stdout.jsonl" \
  2> "${STAGE}/validation/prediction.stderr.txt"
python -m pip freeze --all > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  "${ROOT}/scripts/sei_predict_gse281364.py" \
  "${ROOT}/tests/unit/test_sei_predict_gse281364.py" \
  "${ROOT}/slurm/run_sei_gse281364_smoke.sbatch" \
  > "${STAGE}/source.sha256"
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "predictions/receipt.json").read_text())
if (
    receipt["status"] != "pass_outcome_blind_native_prediction"
    or receipt["elements"] != 1033
    or receipt["outer_locus_sequence_groups"] != 1033
    or receipt["outer_folds"] != 5
    or receipt["native_sequence_classes"] != 40
    or receipt["repeat_max_abs"] != 0.0
    or receipt["outcomes_read"]
    or receipt["reporter_counts_read"]
    or receipt["sealed_outcomes_read"]
    or receipt["model_fitted_or_adapted"]
    or receipt["champion_eligible"]
):
    raise SystemExit("Sei GSE281364 prediction receipt differs")
freeze_tree(root, {
    "artifact_class": "sei_gse281364_outcome_blind_prediction",
    "dataset_id": "gse281364",
    "model_id": "sei",
    "elements": 1033,
    "outer_locus_sequence_groups": 1033,
    "outcomes_read": False,
    "model_fitted_or_adapted": False,
    "champion_eligible": False,
    "status": "passed",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Sei GSE281364 outcome-blind prediction: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
