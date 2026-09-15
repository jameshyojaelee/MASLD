#!/usr/bin/env bash
#SBATCH --job-name=model-work-203
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=16:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export KERAS_BACKEND=tensorflow
export TF_DETERMINISTIC_OPS=1
export TF_CPP_MIN_LOG_LEVEL=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME="${ROOT}/executions/environments/enformer-crested-tensorflow-2.17.0-21065799"
ADMISSION="${ROOT}/executions/enformer-crested-restricted-port-admission-21065579"
INSPECTION="${ROOT}/executions/enformer-crested-safe-keras-inspection-21065759"
CROSSWALK="${ROOT}/executions/enformer-crested-track-crosswalk-21065831"
FORWARD="${ROOT}/executions/enformer-crested-restricted-forward-21066179"
FASTA="${ROOT}/executions/sequence-fasta-21063808"
FIXTURE="${ROOT}/executions/gse281364-enformer-fixture-21070094"
FIXTURE_AUDIT="${ROOT}/executions/gse281364-enformer-fixture-validation-21070139"
PREFLIGHT="${ROOT}/executions/model-work-203-preflight-21070110"
FIXTURE_SHA256=0faf9d48511a76c282bfc412195deebda63503fb5685cd3fb2bfd77c72505231
CACHE="${ROOT}/executions/resumable-model-work-203"
CACHE_ID=restricted-enformer-gse281364-v1
TARGET="${ROOT}/executions/model-work-203-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite restricted Enformer prediction: %s\n' "${TARGET}" >&2
  exit 1
fi
for pair in \
  "3534aee9e54d8fe6dc09c01c63aaead911590b97a08a7a686d56024520885ec6 ${RUNTIME}/ARTIFACTS.json" \
  "e69877c43f47ce7d40962129dae7819a66123ac3f5bef7e0283ac654f77006a5 ${ADMISSION}/ARTIFACTS.json" \
  "155537e853b730d440244581fda83ddff9ba8c13c9dabaa1f57bf841bb69196f ${INSPECTION}/ARTIFACTS.json" \
  "a658a159b03fb05c9a468b4fdbc9c4d745916ad3579b7752fcc4254074dbea10 ${CROSSWALK}/ARTIFACTS.json" \
  "2bb24f274fd58f3011eef9e221e395392594391f57a8a17f8e9ecb7a510771c3 ${FORWARD}/ARTIFACTS.json" \
  "7c202e392c3e415ee7f6bc196b931199c514b2fefc04a19fc82645d81fd8ba39 ${FASTA}/ARTIFACTS.json" \
  "${FIXTURE_SHA256} ${FIXTURE}/ARTIFACTS.json" \
  "20d35712aa9ffbba6708c82731c76b9b4c8b90d214e736646f62c30a4868c198 ${FIXTURE_AUDIT}/ARTIFACTS.json" \
  "edcbedff8b5f0df2d0773855f495b394d5ac24d4eda7f5b76917997ea3bcc810 ${PREFLIGHT}/ARTIFACTS.json"; do
  read -r expected path <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${path}" | sha256sum --check --strict
done
for pair in \
  "4ea75d243ab640e44384b0d79c6d71cdb5e7b6beb2e1a45c297f32d614eb9612 ${ROOT}/scripts/enformer_predict_gse281364.py" \
  "6eca362893d0d64cee64b9a75d399fa1ca70a3a342ff3ae69f02700263ee713f ${ROOT}/scripts/enformer_crested_probe_checkpoint.py" \
  "771c9eeceb5281839e105c47c71a8393b67052633777188b03f336514cc0d68d ${ROOT}/tests/unit/test_enformer_predict_gse281364.py"; do
  read -r expected path <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${path}" | sha256sum --check --strict
done
install -d -m 0750 "${STAGE}/scratch" "${STAGE}/validation" "${CACHE}"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${RUNTIME}" "${ADMISSION}" "${INSPECTION}" "${CROSSWALK}" "${FORWARD}" "${FASTA}" "${FIXTURE}" "${FIXTURE_AUDIT}" "${PREFLIGHT}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for value in sys.argv[1:]:
    verify_frozen_tree(Path(value).resolve(strict=True))
PY
module purge
module load singularity/4.2.2
IMAGE="${RUNTIME}/tensorflow-tensorflow-2.17.0-gpu.sif"
container=(
  singularity exec --nv --cleanenv --containall --net --network none
  --env PYTHONNOUSERSITE=1
  --env PYTHONDONTWRITEBYTECODE=1
  --env PYTHONPATH="${ROOT}"
  --env KERAS_BACKEND=tensorflow
  --env TF_DETERMINISTIC_OPS=1
  --env TF_CPP_MIN_LOG_LEVEL=1
  --bind "${ROOT}:${ROOT}:ro"
  --bind "${STAGE}:${STAGE}:rw"
  --bind "${CACHE}:${CACHE}:rw"
  --pwd "${ROOT}"
  "${IMAGE}"
)
cd "${ROOT}"
"${container[@]}" python -m unittest -v tests.unit.test_enformer_predict_gse281364 \
  > "${STAGE}/validation/tests.stdout.txt" 2> "${STAGE}/validation/tests.stderr.txt"
"${container[@]}" python scripts/enformer_predict_gse281364.py \
  --model-archive "${ADMISSION}/downloads/enformer_human.tar.gz" \
  --source-archive "${ADMISSION}/downloads/crested-source.tar.gz" \
  --fixture-artifacts "${FIXTURE}/ARTIFACTS.json" \
  --manifest "${FIXTURE}/fixture/manifest.tsv" \
  --fasta-artifacts "${FASTA}/ARTIFACTS.json" \
  --fasta "${FASTA}/GRCh38.p14.sequence_model.fa" \
  --fai "${FASTA}/GRCh38.p14.sequence_model.fa.fai" \
  --track-crosswalk-artifacts "${CROSSWALK}/ARTIFACTS.json" \
  --track-crosswalk "${CROSSWALK}/crosswalk/track_crosswalk.tsv" \
  --checkpoint-root "${CACHE}" --cache-id "${CACHE_ID}" \
  --scratch "${STAGE}/scratch" --chunk-size 16 \
  --output "${STAGE}/predictions" \
  > "${STAGE}/validation/prediction.stdout.jsonl" \
  2> "${STAGE}/validation/prediction.stderr.txt"
# The restricted loader materializes exactly three pinned files. Validate their
# complete inventory, type, byte size, and digest before deleting only these
# temporary copies; restricted checkpoint/source bytes must not enter the
# prediction file.
[[ "$(find "${STAGE}/scratch" -mindepth 1 -maxdepth 1 -printf '%f\n' | LC_ALL=C sort)" == \
  $'crested_attention.py\ncrested_layers.py\nenformer_crested_human.keras' ]] || {
  printf 'Enformer scratch inventory differs\n' >&2; exit 1;
}
for path in \
  "${STAGE}/scratch/enformer_crested_human.keras" \
  "${STAGE}/scratch/crested_attention.py" \
  "${STAGE}/scratch/crested_layers.py"; do
  [[ -f "${path}" && ! -L "${path}" ]] || { printf 'Enformer scratch type differs: %s\n' "${path}" >&2; exit 1; }
done
[[ "$(stat -c '%s' "${STAGE}/scratch/enformer_crested_human.keras")" == 985569856 ]] || exit 1
[[ "$(stat -c '%s' "${STAGE}/scratch/crested_attention.py")" == 25834 ]] || exit 1
[[ "$(stat -c '%s' "${STAGE}/scratch/crested_layers.py")" == 23651 ]] || exit 1
for pair in \
  "29dc3835c1d13a6c6bd93315b554075107cc409a9179b11ee0e00180f67cef56 ${STAGE}/scratch/enformer_crested_human.keras" \
  "6c195deb31303d5c1664341308f84625c3d862038fafc8725c781ce33d1c0dac ${STAGE}/scratch/crested_attention.py" \
  "44d4befd0e7faca5bf4eea57e94c296147d51818b805d6a8a8a2189ef2f4b457 ${STAGE}/scratch/crested_layers.py"; do
  read -r expected path <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${path}" | sha256sum --check --strict
done
rm -- \
  "${STAGE}/scratch/enformer_crested_human.keras" \
  "${STAGE}/scratch/crested_attention.py" \
  "${STAGE}/scratch/crested_layers.py"
rmdir "${STAGE}/scratch"
"${container[@]}" python -m pip freeze --all > "${STAGE}/validation/effective-pip-freeze.txt"
"${container[@]}" nvidia-smi -q > "${STAGE}/validation/nvidia-smi.txt"
module list 2> "${STAGE}/validation/module-list.txt"
sha256sum scripts/enformer_predict_gse281364.py \
  scripts/enformer_crested_probe_checkpoint.py \
  tests/unit/test_enformer_predict_gse281364.py \
  slurm/run_model_work_203.sbatch > "${STAGE}/source.sha256"
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
    receipt["status"] != "pass_outcome_blind_restricted_prediction"
    or receipt["elements"] != 1033
    or receipt["outer_locus_sequence_groups"] != 1033
    or receipt["outer_folds"] != 5
    or receipt["tracks"] != 5313
    or receipt["accessibility_tracks"] != 684
    or receipt["deterministic_repeat_max_abs"] != 0.0
    or receipt["outcomes_read"]
    or receipt["reporter_counts_read"]
    or receipt["sealed_outcomes_read"]
    or receipt["model_fitted_or_adapted"]
    or receipt["native_sonnet_parity_established"]
    or receipt["champion_eligible"]
):
    raise SystemExit("restricted Enformer prediction receipt differs")
freeze_tree(root, {
    "artifact_class": "enformer_gse281364_outcome_blind_prediction",
    "dataset_id": "gse281364",
    "model_id": "enformer_crested_restricted_port",
    "elements": 1033,
    "outer_locus_sequence_groups": 1033,
    "tracks": 5313,
    "outcomes_read": False,
    "model_fitted_or_adapted": False,
    "native_sonnet_parity_established": False,
    "champion_eligible": False,
    "status": "passed",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Restricted Enformer GSE281364 prediction: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
