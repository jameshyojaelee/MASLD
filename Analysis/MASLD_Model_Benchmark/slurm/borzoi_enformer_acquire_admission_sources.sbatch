#!/usr/bin/env bash
#SBATCH --job-name=masld-borzoi-enformer-admit
#SBATCH --partition=io
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
CONFIG="${ROOT}/config/borzoi_enformer_admission.json"
TARGET="${ROOT}/executions/borzoi-enformer-admission-sources-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite Borzoi/Enformer admission source artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/sources" "${STAGE}/validation"
curl_flags=(--fail --location --retry 6 --retry-all-errors --connect-timeout 30)

BORZOI_REV=9736924e6c861f0d71978cc3012ff8d76d6fa91f
ENFORMER_REV=f5de0ede8430809180254ee957abf36ed62579ef
curl "${curl_flags[@]}" --output "${STAGE}/sources/borzoi_license" \
  "https://raw.githubusercontent.com/calico/borzoi/${BORZOI_REV}/LICENSE"
curl "${curl_flags[@]}" --output "${STAGE}/sources/borzoi_readme" \
  "https://raw.githubusercontent.com/calico/borzoi/${BORZOI_REV}/README.md"
curl "${curl_flags[@]}" --output "${STAGE}/sources/borzoi_params_pred" \
  "https://raw.githubusercontent.com/calico/borzoi/${BORZOI_REV}/examples/params_pred.json"
curl "${curl_flags[@]}" --output "${STAGE}/sources/borzoi_targets_human.txt.gz" \
  "https://raw.githubusercontent.com/calico/borzoi/${BORZOI_REV}/data/targets_human.txt.gz"
curl "${curl_flags[@]}" --output "${STAGE}/sources/borzoi_checkpoint_params.json" \
  'https://storage.googleapis.com/download/storage/v1/b/seqnn-share/o/borzoi%2Fparams.json?alt=media&generation=1699567843281710'

curl "${curl_flags[@]}" --output "${STAGE}/sources/enformer_license" \
  "https://raw.githubusercontent.com/google-deepmind/deepmind-research/${ENFORMER_REV}/LICENSE"
curl "${curl_flags[@]}" --output "${STAGE}/sources/enformer_readme" \
  "https://raw.githubusercontent.com/google-deepmind/deepmind-research/${ENFORMER_REV}/enformer/README.md"
curl "${curl_flags[@]}" --output "${STAGE}/sources/enformer_py" \
  "https://raw.githubusercontent.com/google-deepmind/deepmind-research/${ENFORMER_REV}/enformer/enformer.py"
curl "${curl_flags[@]}" --output "${STAGE}/sources/enformer_attention" \
  "https://raw.githubusercontent.com/google-deepmind/deepmind-research/${ENFORMER_REV}/enformer/attention_module.py"
curl "${curl_flags[@]}" --output "${STAGE}/sources/enformer_requirements" \
  "https://raw.githubusercontent.com/google-deepmind/deepmind-research/${ENFORMER_REV}/enformer/requirements.txt"
curl "${curl_flags[@]}" --output "${STAGE}/sources/enformer_targets_human.txt" \
  'https://raw.githubusercontent.com/calico/basenji/0.5/manuscripts/cross2020/targets_human.txt'

for index in 0 1 2 3; do
  case "${index}" in
    0) object='borzoi%2Ff0%2Fmodel0_best.h5'; generation=1691539358988138 ;;
    1) object='borzoi%2Ff1%2Fmodel0_best.h5'; generation=1691539419241923 ;;
    2) object='borzoi%2Ff2%2Fmodel0_best.h5'; generation=1691539441571542 ;;
    3) object='borzoi%2Ff3%2Fmodel0_best.h5'; generation=1691539462862235 ;;
  esac
  curl "${curl_flags[@]}" \
    --output "${STAGE}/sources/borzoi_checkpoint_${index}.metadata.json" \
    "https://storage.googleapis.com/storage/v1/b/seqnn-share/o/${object}?generation=${generation}"
done
enformer_objects=(
  'models%2Fenformer%2Fsonnet_weights%2Fcheckpoint'
  'models%2Fenformer%2Fsonnet_weights%2Fenformer-fine-tuned-human-1.data-00000-of-00001'
  'models%2Fenformer%2Fsonnet_weights%2Fenformer-fine-tuned-human-1.index'
)
enformer_generations=(1653413781220972 1653414053766013 1653414056026643)
for index in 0 1 2; do
  curl "${curl_flags[@]}" \
    --output "${STAGE}/sources/enformer_checkpoint_${index}.metadata.json" \
    "https://storage.googleapis.com/storage/v1/b/dm-enformer/o/${enformer_objects[$index]}?generation=${enformer_generations[$index]}"
done

gzip -t "${STAGE}/sources/borzoi_targets_human.txt.gz"
module load python/3.11.5-GCCcore-13.2.0
python "${ROOT}/scripts/borzoi_enformer_validate_sources.py" \
  --config "${CONFIG}" \
  --source-dir "${STAGE}/sources" \
  --project-root "${ROOT}" \
  --output "${STAGE}/review" \
  > "${STAGE}/validation/source-validation.stdout.txt" \
  2> "${STAGE}/validation/source-validation.stderr.txt"
curl --version > "${STAGE}/validation/curl.version.txt"
python --version > "${STAGE}/validation/python.version.txt" 2>&1
python -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
sha256sum \
  "${CONFIG}" \
  "${ROOT}/scripts/borzoi_enformer_validate_sources.py" \
  "${ROOT}/slurm/borzoi_enformer_acquire_admission_sources.sbatch" \
  > "${STAGE}/source.sha256"

PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "review" / "source_admission_receipt.json").read_text())
freeze_tree(
    root,
    {
        "artifact_class": "borzoi_enformer_admission_sources",
        "models": ["borzoi_ensemble", "enformer"],
        "checkpoint_bytes_downloaded": receipt["checkpoint_bytes_downloaded"],
        "model_forward_executed": receipt["model_forward_executed"],
        "terminal_disposition": receipt["terminal_disposition"],
        "independent_source_license_exposure_review": True,
        "status": "passed",
    },
)
PY
mv "${STAGE}" "${TARGET}"
printf 'Borzoi/Enformer source admission artifact: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
