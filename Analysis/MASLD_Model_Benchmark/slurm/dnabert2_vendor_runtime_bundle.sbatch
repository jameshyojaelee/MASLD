#!/usr/bin/env bash
#SBATCH --job-name=masld-dnabert2-code
#SBATCH --partition=io
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
SOURCE="${ROOT}/executions/dna-lm-admission-sources-21064540"
SOURCE_SHA=4e1703b03ff28924904053e842b89079d510a80d2a9ad082de5369adaeea75c8
REVISION=7bce263b15377fc15361f52cfab88f8b586abda0
TARGET="${ROOT}/executions/dnabert2-runtime-bundle-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite DNABERT-2 runtime bundle: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/bundle" "${STAGE}/validation"
printf '%s  %s\n' "${SOURCE_SHA}" "${SOURCE}/ARTIFACTS.json" \
  | sha256sum --check --strict

BASE_URL="https://huggingface.co/zhihan1996/DNABERT-2-117M/resolve/${REVISION}"
for filename in \
  LICENSE \
  bert_layers.py \
  bert_padding.py \
  config.json \
  configuration_bert.py \
  flash_attn_triton.py \
  tokenizer.json \
  tokenizer_config.json; do
  curl --fail --location --retry 6 --retry-all-errors \
    --output "${STAGE}/bundle/${filename}" "${BASE_URL}/${filename}"
done

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${SOURCE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

verify_frozen_tree(Path(sys.argv[1]))
PY
python "${ROOT}/scripts/dnabert2_validate_runtime_bundle.py" \
  --source "${SOURCE}" \
  --bundle "${STAGE}/bundle" \
  --output "${STAGE}/review" \
  > "${STAGE}/validation/validator.stdout.txt" \
  2> "${STAGE}/validation/validator.stderr.txt"
python -m pip freeze > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  "${ROOT}/scripts/dnabert2_validate_runtime_bundle.py" \
  "${ROOT}/slurm/dnabert2_vendor_runtime_bundle.sbatch" \
  "${SOURCE}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"

PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "review" / "runtime_bundle_receipt.json").read_text())
if (
    receipt["status"] != "pass"
    or receipt["custom_code_executed"]
    or receipt["checkpoint_loaded"]
    or receipt["model_forward_executed"]
):
    raise SystemExit("DNABERT-2 runtime bundle receipt differs")
freeze_tree(
    root,
    {
        "artifact_class": "dnabert2_revision_pinned_runtime_bundle",
        "upstream_source_artifacts_sha256": "4e1703b03ff28924904053e842b89079d510a80d2a9ad082de5369adaeea75c8",
        "checkpoint_revision": "7bce263b15377fc15361f52cfab88f8b586abda0",
        "license": "apache-2.0",
        "custom_code_executed": False,
    },
)
PY

mv "${STAGE}" "${TARGET}"
printf 'DNABERT-2 runtime bundle: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
