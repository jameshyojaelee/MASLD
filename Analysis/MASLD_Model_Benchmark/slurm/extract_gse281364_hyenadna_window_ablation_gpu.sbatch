#!/usr/bin/env bash
#SBATCH --job-name=hyenadna-window-ablation-extract
#SBATCH --partition=gpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --gres=gpu:l40s:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export CUBLAS_WORKSPACE_CONFIG=:4096:8
export TOKENIZERS_PARALLELISM=false
export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
ABLATION="${ROOT}/executions/gse281364-hyenadna-window-ablation-20260905T192600Z"
PYTHON=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/.agent/scanpy_env/bin/python
RUNTIME="${ROOT}/executions/environments/dna-lm-transformers-4.38.1-21064683"
CODE="${ROOT}/executions/dna-lm-custom-code-21064626"
CONFIG="${ROOT}/config/gse281364_dna_lm_native_contract.json"
OUTPUT="${ABLATION}/embeddings"

install -d -m 0750 "${ABLATION}/validation"
if [[ -e "${OUTPUT}" ]]; then
  printf 'Refusing to overwrite window-ablation embeddings: %s\n' "${OUTPUT}" >&2
  exit 1
fi

# The checkpoint, custom code, and 4,096 bp fixture must be the frozen ones the
# published result was extracted from.
printf '%s  %s\n' \
  "70906faae158909a8972f674f680943d8a51467032ae3266a6a12f2152beb36c" "${CODE}/ARTIFACTS.json" \
  "6d21df999835624905e08573038dc72e893a5b25a8ca8f4ec1acea8f61d91c77" "${ROOT}/executions/dna-lm-open-safetensors-21064556/ARTIFACTS.json" \
  "74a3a5c36bd1416daf3fe4b63fe317cfa7e2d008b96e2503d4ca43770b4357de" "${ROOT}/executions/gse281364-dna-lm-common-fixture-21069076/ARTIFACTS.json" \
  | sha256sum --check --strict

export PYTHONPATH="${ROOT}:${RUNTIME}/site-packages:${CODE}/code"
nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv,noheader \
  > "${ABLATION}/validation/nvidia-smi.csv"

"${PYTHON}" "${ROOT}/scripts/extract_gse281364_hyenadna_window_ablation.py" \
  --root "${ROOT}" \
  --config "${CONFIG}" \
  --output "${OUTPUT}" \
  --seed 20260824 \
  --causality-probe-elements 8 \
  > "${ABLATION}/validation/extraction.stdout.json" \
  2> "${ABLATION}/validation/extraction.stderr.txt"

"${PYTHON}" -m pip freeze > "${ABLATION}/validation/extraction.pip-freeze.txt"
"${PYTHON}" --version > "${ABLATION}/validation/extraction.python.version.txt" 2>&1
sha256sum \
  "${ROOT}/scripts/extract_gse281364_hyenadna_window_ablation.py" \
  "${ROOT}/scripts/gse281364_dna_lm_native_contract.py" \
  "${CONFIG}" \
  "${ABLATION}/prespec/PRESPECIFICATION.md" \
  > "${ABLATION}/validation/extraction.source.sha256"

printf 'Completed HyenaDNA window-ablation extraction: %s\n' "${OUTPUT}"
