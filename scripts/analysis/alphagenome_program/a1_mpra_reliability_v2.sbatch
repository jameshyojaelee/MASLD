#!/bin/bash
#SBATCH --job-name=agp_a1-mpra-reliability-v2
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_a1_reliability_v2_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_a1_reliability_v2_%j.err

set -euo pipefail
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export A1_N_BOOT=10000
export A1_N_PERM=10000

PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SRC="$PROJ/scripts/analysis/alphagenome_program"
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
OUT="$1"

echo "host $(hostname)  start $(date -u +%Y-%m-%dT%H:%M:%SZ)  jobid ${SLURM_JOB_ID}"
echo "A1_N_BOOT=${A1_N_BOOT}  A1_N_PERM=${A1_N_PERM}"
mkdir -p "$OUT"
cp "$SRC/a1_mpra_reliability.py"    "$OUT/a1_mpra_reliability.py.run"
cp "$SRC/a1_mpra_reliability_v2.py" "$OUT/a1_mpra_reliability_v2.py.run"
cp "$SRC/a1_mpra_reliability_v2.sbatch" "$OUT/a1_mpra_reliability_v2.sbatch.run"
sha256sum "$SRC/a1_mpra_reliability.py" "$SRC/a1_mpra_reliability_v2.py" \
          "$SRC/a1_mpra_reliability_v2.sbatch" > "$OUT/script.sha256"

echo "=== stage 1: full wave-1 analysis, 10000 draws / 10000 permutations ==="
"$PY" "$SRC/a1_mpra_reliability.py" "$OUT"

echo "=== stage 2: v2 addenda (the eight corrections) ==="
"$PY" "$SRC/a1_mpra_reliability_v2.py" "$OUT"

echo "done $(date -u +%Y-%m-%dT%H:%M:%SZ)"
