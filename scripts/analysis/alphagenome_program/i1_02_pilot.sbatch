#!/bin/bash
#SBATCH --job-name=agp_ad1_pilot
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_ad1_pilot_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_ad1_pilot_%j.err

# AD.1 cost-and-length pilot. The prespecification requires compute to be measured before it is
# spent, and the archived Atlas never recorded the sequence length it scored at, so the same 64
# leads are scored at all five SDK lengths and the length whose scores match the archive is the one
# the full gate runs at. L40S only: scores are not reproducible across GPU architectures.

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

mkdir -p "$RUN/raw" "$RUN/env" "$RUN/logs" "$RUN/tables"

echo "=== node $(hostname) job $SLURM_JOB_ID $(date -u +%Y-%m-%dT%H:%M:%SZ) ==="
nvidia-smi --query-gpu=name,memory.total,driver_version,compute_cap --format=csv \
  | tee "$RUN/env/gpu_pilot.csv"

$MM run -n $ENV pip freeze > "$RUN/env/pip_freeze.txt" 2>/dev/null
$MM run -n $ENV python -c "import sys; print(sys.version)" > "$RUN/env/python_version.txt" 2>&1
$MM run -n $ENV python - <<'PY' > "$RUN/env/versions.txt" 2>&1
import sys, platform, jax, jaxlib, numpy, scipy, pandas
import alphagenome, alphagenome_research
print("python", sys.version.replace("\n", " "))
print("platform", platform.platform())
print("jax", jax.__version__); print("jaxlib", jaxlib.__version__)
print("numpy", numpy.__version__); print("scipy", scipy.__version__)
print("pandas", pandas.__version__)
print("alphagenome", getattr(alphagenome, "__version__", "n/a"))
print("alphagenome_research", getattr(alphagenome_research, "__version__", "n/a"))
print("jax devices", jax.devices()); print("backend", jax.default_backend())
PY
cat "$RUN/env/versions.txt"

# ---- length sweep: the same 64 leads, one per 1-Mb block, at every SDK length
for L in 2048 16384 131072 524288 1048576; do
  echo "=== sweep length $L ==="
  $MM run -n $ENV python "$SRC/i1_01_score_local.py" \
      --variants "$RUN/inputs/sweep_variants.tsv" \
      --length "$L" \
      --out "$RUN/raw/sweep_${L}.tsv"
  echo "  exit=$?"
done

# ---- cost pilot at 1 Mb on the block-spread pilot subset, for the full-run time estimate
echo "=== cost pilot 1048576 ==="
$MM run -n $ENV python "$SRC/i1_01_score_local.py" \
    --variants "$RUN/inputs/pilot_variants.tsv" \
    --length 1048576 \
    --out "$RUN/raw/pilot_1048576.tsv"
echo "  exit=$?"

# ---- read the sweep against the archive and pick the length
echo "=== sweep analysis ==="
$MM run -n $ENV python "$SRC/i1_03_pilot_analysis.py" --run "$RUN"
echo "  exit=$?"

echo "=== done $(date -u +%Y-%m-%dT%H:%M:%SZ) ==="
