#!/bin/bash
#SBATCH --job-name=agp_ad2_mpra_embed
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_ad2_mpra_embed_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_ad2_mpra_embed_%j.err

set -uo pipefail

export PYTHONNOUSERSITE=1
export PYTHONDONTWRITEBYTECODE=1
export MAMBA_ROOT_PREFIX=/gpfs/commons/home/jameslee/micromamba
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export TF_CPP_MIN_LOG_LEVEL=3
export HF_HUB_OFFLINE=1

MM=/gpfs/commons/home/jameslee/.local/bin/micromamba
ENV=alphagenome_local
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SRC=$PROJ/scripts/analysis/alphagenome_program
OUT=$PROJ/GWAS/finemapping/results/alphagenome_program/ad-arm2-mpra-20260915T114407Z

mkdir -p "$OUT/embeddings" "$OUT/env"

echo "=== node $(hostname) job ${SLURM_JOB_ID} $(date -u +%Y-%m-%dT%H:%M:%SZ) ==="
nvidia-smi | tee "$OUT/env/nvidia_smi.txt"
nvidia-smi --query-gpu=name,memory.total,driver_version,compute_cap --format=csv \
  | tee "$OUT/env/gpu.csv"

# This arm scores on L40S only. Amendment 01 section 3: the same weights and sequence gave a 0.0388
# channel gap across GPU architectures, the size of a typical measured indel effect.
GPUNAME=$(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)
case "$GPUNAME" in
  *L40S*) CARD=l40s ;;
  *)      echo "REFUSING: this arm requires an L40S, the node offers '$GPUNAME'" >&2; exit 1 ;;
esac
echo "card=$CARD name=$GPUNAME"

$MM run -n $ENV pip freeze > "$OUT/env/pip_freeze_gpu.txt" 2>/dev/null
$MM run -n $ENV python - > "$OUT/env/versions_gpu.txt" 2>&1 <<'PY'
import sys, platform, jax, jaxlib, numpy, scipy
import alphagenome, alphagenome_research
print("python", sys.version.replace("\n", " "))
print("platform", platform.platform())
print("jax", jax.__version__, "jaxlib", jaxlib.__version__)
print("numpy", numpy.__version__, "scipy", scipy.__version__)
print("alphagenome", getattr(alphagenome, "__version__", "n/a"))
print("alphagenome_research", getattr(alphagenome_research, "__version__", "n/a"))
print("devices", jax.local_devices())
print("device_kind", [d.device_kind for d in jax.local_devices()])
PY
cat "$OUT/env/versions_gpu.txt"

RC=0
for L in 2048 16384; do
  echo "=== extract ${L} bp ==="
  $MM run -n $ENV python "$SRC/i1_extract_mpra_embeddings.py" \
      --length "$L" \
      --card "$CARD" \
      --out "$OUT/embeddings/alphagenome_frozen_trunk_${L}bp.npz" \
      --receipt "$OUT/embeddings/extraction_receipt_${L}bp.json"
  S=$?
  echo "  exit=${S}"
  [ "$S" -ne 0 ] && RC=$S
done

echo "=== done $(date -u +%Y-%m-%dT%H:%M:%SZ) rc=${RC} ==="
exit $RC
