#!/bin/bash
#SBATCH --job-name=agp_ad1_gate
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=14:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_ad1_gate_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_ad1_gate_%j.err

# AD.1 full gate: all 3,845 matched Currin leads, REF and ALT, at the length the pilot chose.
# L40S only, because the runtime probe measured that scores are bit-identical on one card and not
# reproducible across GPU architectures by as much as a typical measured effect.
# --resume makes a wall-time kill recoverable: the score table is appended and flushed every 100
# variants, and a restart skips the keys already on disk for this card.

set -uo pipefail

export PYTHONNOUSERSITE=1
export MAMBA_ROOT_PREFIX=/gpfs/commons/home/jameslee/micromamba
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export TF_CPP_MIN_LOG_LEVEL=3
export HF_HUB_OFFLINE=1

MM=/gpfs/commons/home/jameslee/.local/bin/micromamba
ENV=alphagenome_local
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SRC=$PROJ/scripts/analysis/alphagenome_program
RUN=$PROJ/GWAS/finemapping/results/alphagenome_program/ad1-gate-20260915T114615Z
LENGTH=${AD1_LENGTH:?set AD1_LENGTH from the pilot verdict before submitting}

echo "=== node $(hostname) job $SLURM_JOB_ID $(date -u +%Y-%m-%dT%H:%M:%SZ) length $LENGTH ==="
nvidia-smi --query-gpu=name,memory.total,driver_version,compute_cap --format=csv \
  | tee "$RUN/env/gpu_gate.csv"

$MM run -n $ENV python "$SRC/i1_01_score_local.py" \
    --variants "$RUN/inputs/gate_variants.tsv" \
    --length "$LENGTH" \
    --out "$RUN/raw/gate_${LENGTH}.tsv" \
    --resume
echo "  score exit=$?"

echo "=== gate statistics ==="
$MM run -n $ENV python "$SRC/i1_04_gate_stats.py" \
    --run "$RUN" --scores "$RUN/raw/gate_${LENGTH}.tsv"
echo "  exit=$?"

echo "=== done $(date -u +%Y-%m-%dT%H:%M:%SZ) ==="
