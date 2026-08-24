#!/usr/bin/env bash
#SBATCH --job-name=masld-dna-lm-code
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
TARGET="${ROOT}/executions/dna-lm-custom-code-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite DNA-LM custom-code artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/code/hyenadna" "${STAGE}/code/caduceus" "${STAGE}/validation"
printf '%s  %s\n' 4e1703b03ff28924904053e842b89079d510a80d2a9ad082de5369adaeea75c8 \
  "${SOURCE}/ARTIFACTS.json" | sha256sum --check --strict

curl_flags=(--fail --location --retry 6 --retry-all-errors --connect-timeout 30)
HYENA_REPO="https://huggingface.co/LongSafari/hyenadna-large-1m-seqlen-hf/resolve/0a629abf9c7f85b4ec9aa6a1aefa3adcf1907446"
for filename in config.json configuration_hyena.py modeling_hyena.py tokenization_hyena.py; do
  curl "${curl_flags[@]}" --output "${STAGE}/code/hyenadna/${filename}" "${HYENA_REPO}/${filename}"
done
CADUCEUS_REPO="https://huggingface.co/kuleshov-group/caduceus-ps_seqlen-131k_d_model-256_n_layer-16/resolve/d89eeb853136ea64da7feb3d0c8e909771b17ae6"
for filename in config.json configuration_caduceus.py modeling_caduceus.py modeling_rcps.py tokenization_caduceus.py; do
  curl "${curl_flags[@]}" --output "${STAGE}/code/caduceus/${filename}" "${CADUCEUS_REPO}/${filename}"
done

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${SOURCE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

verify_frozen_tree(Path(sys.argv[1]))
print("DNA-LM source artifact verified")
PY
python "${ROOT}/scripts/dna_lm_validate_custom_code.py" \
  --source-artifact "${SOURCE}" \
  --code-root "${STAGE}/code" \
  --output "${STAGE}/review" \
  > "${STAGE}/validation/custom-code.stdout.txt" \
  2> "${STAGE}/validation/custom-code.stderr.txt"
curl --version > "${STAGE}/validation/curl.version.txt"
python --version > "${STAGE}/validation/python.version.txt" 2>&1
python -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
sha256sum \
  "${ROOT}/scripts/dna_lm_validate_custom_code.py" \
  "${ROOT}/slurm/dna_lm_vendor_custom_code.sbatch" \
  "${SOURCE}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"

PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "review" / "custom_code_admission_receipt.json").read_text())
if (
    receipt["status"] != "pass"
    or receipt["custom_code_executed"]
    or receipt["checkpoint_tensor_data_loaded"]
    or receipt["model_forward_executed"]
    or receipt["runtime_network_allowed"]
    or receipt["observed_outcomes_loaded"]
):
    raise SystemExit("custom-code admission receipt differs")
freeze_tree(
    root,
    {
        "artifact_class": "dna_lm_custom_code",
        "models": ["hyenadna", "caduceus"],
        "custom_code_executed": False,
        "checkpoint_tensor_data_loaded": False,
        "model_forward_executed": False,
        "status": "passed",
    },
)
PY
mv "${STAGE}" "${TARGET}"
printf 'DNA-LM custom-code artifact: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
