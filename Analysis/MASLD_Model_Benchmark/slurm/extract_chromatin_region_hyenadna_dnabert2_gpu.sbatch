#!/usr/bin/env bash
#SBATCH --job-name=chromatin-region-emb-hyena-dnabert2
#SBATCH --partition=gpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --gres=gpu:l40s:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export TOKENIZERS_PARALLELISM=false
export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
SEQLANE="${ROOT}/executions/rerun-floor-20260905T172141Z/chromatin/seqlane"
PAIRS="${SEQLANE}/pairs_w1024.npz"
# Python 3.10 is required: the transformers 4.38.1 shim ships a `tokenizers`
# extension compiled for cp310, which is what breaks under a 3.11 interpreter.
PYTHON=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/.agent/scanpy_env/bin/python
RUNTIME="${ROOT}/executions/environments/dna-lm-transformers-4.38.1-21064683"
EINOPS="${ROOT}/executions/environments/dnabert2-einops-0.6.1-21064904"
CODE="${ROOT}/executions/dna-lm-custom-code-21064626"
BUNDLE="${ROOT}/executions/dnabert2-runtime-bundle-21064834"

export PYTHONPATH="${ROOT}:${EINOPS}/site-packages:${RUNTIME}/site-packages:${CODE}/code:${BUNDLE}"
nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv,noheader

for model in hyenadna dnabert2; do
  for mode in pretrained random; do
    out="${SEQLANE}/emb_${model}_${mode}.npz"
    if [[ -e "${out}" ]]; then
      printf 'Skipping existing %s\n' "${out}"
      continue
    fi
    "${PYTHON}" "${ROOT}/scripts/extract_chromatin_region_dna_lm_embeddings.py" \
      --root "${ROOT}" \
      --pairs "${PAIRS}" \
      --model "${model}" \
      --mode "${mode}" \
      --output "${out}" \
      --seed 20260905 \
      --batch 32
  done
done

"${PYTHON}" -m pip freeze > "${SEQLANE}/emb_extraction.pip-freeze.txt"
printf 'Completed HyenaDNA and DNABERT-2 region embeddings\n'
