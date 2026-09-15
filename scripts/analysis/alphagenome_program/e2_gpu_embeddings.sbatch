#!/usr/bin/env bash
#SBATCH --job-name=agp_e2-transfer-gpu
#SBATCH --partition=gpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --gres=gpu:l40s:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_e2-gpu_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_e2-gpu_%j.err

set -euo pipefail
export PYTHONNOUSERSITE=1
export PYTHONDONTWRITEBYTECODE=1
export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1
export CUBLAS_WORKSPACE_CONFIG=:4096:8
export TOKENIZERS_PARALLELISM=false

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
BENCH="${ROOT}/Analysis/MASLD_Model_Benchmark"
OUT="${ROOT}/GWAS/finemapping/results/alphagenome_program/e2-transfer-20260914T203752Z"
SCRIPT="${ROOT}/scripts/analysis/alphagenome_program/c2_02_extract_embeddings.py"
RUNTIME="${BENCH}/executions/environments/dna-lm-transformers-4.38.1-21064683/site-packages"

nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv,noheader

# Same extractor, same checkpoint, same 4,096-bp geometry as the published reporter arm
# w4096_proportional and as C2; only the window table changes.
export PYTHONPATH="${RUNTIME}"
"${ROOT}/.agent/scanpy_env/bin/python" "${SCRIPT}" \
  --tokens "${OUT}/inputs/window_tokens.npy" \
  --meta "${OUT}/inputs/window_meta.npz" \
  --model hyenadna \
  --output "${OUT}/raw/hyenadna_allele_embeddings.npz" \
  --seed 20260914 --batch 16
"${ROOT}/.agent/scanpy_env/bin/python" -m pip freeze \
  > "${OUT}/validation/pip_freeze_hyenadna_runtime.txt"

printf 'agp e2 embeddings complete\n'
