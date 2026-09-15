#!/usr/bin/env bash
#SBATCH --job-name=chromatin-region-emb-caduceus
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
export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
SEQLANE="${ROOT}/executions/rerun-floor-20260905T172141Z/chromatin/seqlane"
PAIRS="${SEQLANE}/pairs_w1024.npz"
RUNTIME="${ROOT}/executions/environments/caduceus-native-torch2.2.0-21065466"
CODE="${ROOT}/executions/dna-lm-custom-code-21064626"
# Caduceus needs mamba_ssm, which lives only in the packed native runtime.
RUNTIME_TMP="${TMPDIR:-/tmp}/caduceus-native-${SLURM_JOB_ID}"
install -d -m 0750 "${RUNTIME_TMP}"
trap 'rm -rf "${RUNTIME_TMP}"' EXIT
tar -xzf "${RUNTIME}/payload/caduceus-native-env.tar.gz" -C "${RUNTIME_TMP}"

export PYTHONPATH="${ROOT}:${CODE}/code"
nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv,noheader
"${RUNTIME_TMP}/bin/python" -c "import torch, mamba_ssm; print('torch', torch.__version__, 'cuda', torch.cuda.is_available())"

for mode in pretrained random; do
  out="${SEQLANE}/emb_caduceus_${mode}.npz"
  if [[ -e "${out}" ]]; then
    printf 'Skipping existing %s\n' "${out}"
    continue
  fi
  "${RUNTIME_TMP}/bin/python" "${ROOT}/scripts/extract_chromatin_region_dna_lm_embeddings.py" \
    --root "${ROOT}" \
    --pairs "${PAIRS}" \
    --model caduceus \
    --mode "${mode}" \
    --output "${out}" \
    --seed 20260905 \
    --batch 16
done

"${RUNTIME_TMP}/bin/python" -m pip freeze > "${SEQLANE}/emb_caduceus.pip-freeze.txt"
printf 'Completed Caduceus region embeddings\n'
