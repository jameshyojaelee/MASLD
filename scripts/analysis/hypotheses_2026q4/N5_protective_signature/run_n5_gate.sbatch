#!/bin/bash
#SBATCH --job-name=N5_gate
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/N5-20261009T201448Z/logs/gate-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/N5-20261009T201448Z/logs/gate-%j.err
# N5 positive-control gate: presence check of HSD17B13 and PNPLA3 perturbations. Metadata only.
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
RUN=$ROOT/Analysis/Hypotheses/N5-20261009T201448Z
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
OUT=$RUN/results/gate-$(date -u +%Y%m%dT%H%M%SZ)
"$PY" -I "$ROOT/scripts/analysis/hypotheses_2026q4/N5_protective_signature/n5_gate_presence_check.py" "$OUT"
"$PY" -I --version > "$OUT/python_version.txt" 2>&1
"$PY" -I -m pip freeze > "$OUT/pip_freeze.txt" 2>&1
sha256sum "$ROOT/scripts/analysis/hypotheses_2026q4/N5_protective_signature/n5_gate_presence_check.py" > "$OUT/script_sha256.txt"
echo "$OUT" > "$RUN/results/LATEST_GATE.txt"
