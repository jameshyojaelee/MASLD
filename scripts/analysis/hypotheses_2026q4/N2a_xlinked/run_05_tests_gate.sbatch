#!/bin/bash
#SBATCH --job-name=N2a_gate
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/N2a-20261009T202200Z/logs/slurm-%x-%j.out
# Synthetic cross-fit test, then the power gate (depth, het calls and autosomal deltas only).
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/N2a-20261009T202200Z
S=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/hypotheses_2026q4/N2a_xlinked
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
sha256sum "$S/05_gate.py" "$S/n2a_lib.py" "$S/test_crossfit.py" >> "$ROOT/logs/script_sha256_05.txt"
"$PY" -I "$S/test_crossfit.py"
OUT="$ROOT/results/gate-$(date -u +%Y%m%dT%H%M%SZ)"
mkdir -p "$ROOT/results"
"$PY" -I "$S/05_gate.py" "$ROOT" "$OUT"
