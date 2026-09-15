#!/usr/bin/env bash
#SBATCH --job-name=agp_c2-endogenous-head-gpu
#SBATCH --partition=gpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --gres=gpu:l40s:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_c2-gpu_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_c2-gpu_%j.err

set -euo pipefail
export PYTHONNOUSERSITE=1
export PYTHONDONTWRITEBYTECODE=1
export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1
export CUBLAS_WORKSPACE_CONFIG=:4096:8
export TOKENIZERS_PARALLELISM=false

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
BENCH="${ROOT}/Analysis/MASLD_Model_Benchmark"
OUT="${ROOT}/GWAS/finemapping/results/alphagenome_program/c2-endogenous-head-20260914T194000Z"
SCRIPT="${ROOT}/scripts/analysis/alphagenome_program/c2_02_extract_embeddings.py"
TOKENS="${OUT}/inputs/window_tokens.npy"
META="${OUT}/inputs/window_meta.npz"
RUNTIME="${BENCH}/executions/environments/dna-lm-transformers-4.38.1-21064683/site-packages"
CADUCEUS_ENV="${BENCH}/executions/environments/caduceus-native-torch2.2.0-21065466"

nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv,noheader

# HyenaDNA: the reporter's runtime (torch 2.3.1 + transformers 4.38.1 overlay), no mamba needed.
export PYTHONPATH="${RUNTIME}"
"${ROOT}/.agent/scanpy_env/bin/python" "${SCRIPT}" \
  --tokens "${TOKENS}" --meta "${META}" \
  --model hyenadna --output "${OUT}/raw/hyenadna_allele_embeddings.npz" \
  --seed 20260914 --batch 16
"${ROOT}/.agent/scanpy_env/bin/python" -m pip freeze > "${OUT}/validation/pip_freeze_hyenadna_runtime.txt"

# Caduceus needs mamba_ssm, which lives only in the packed native runtime.
RUNTIME_TMP="${TMPDIR:-/tmp}/agp-c2-caduceus-${SLURM_JOB_ID}"
install -d -m 0750 "${RUNTIME_TMP}"
trap 'rm -rf "${RUNTIME_TMP}"' EXIT
tar -xzf "${CADUCEUS_ENV}/payload/caduceus-native-env.tar.gz" -C "${RUNTIME_TMP}"
unset PYTHONPATH
"${RUNTIME_TMP}/bin/python" -c "import torch, mamba_ssm; print('torch', torch.__version__, 'cuda', torch.cuda.is_available())"
"${RUNTIME_TMP}/bin/python" "${SCRIPT}" \
  --tokens "${TOKENS}" --meta "${META}" \
  --model caduceus --output "${OUT}/raw/caduceus_allele_embeddings.npz" \
  --seed 20260914 --batch 8
"${RUNTIME_TMP}/bin/python" -m pip freeze > "${OUT}/validation/pip_freeze_caduceus_runtime.txt"

printf 'agp c2 embeddings complete\n'
