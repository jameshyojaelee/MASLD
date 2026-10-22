#!/bin/bash
#SBATCH --job-name=H30_gate
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=90:00:00
#SBATCH --output=%x-%j.out
# H30 detectability gate: depth and background only. Usage: sbatch run_gate.sbatch H30_DIR
set -euo pipefail
H30="${1:?H30 result dir required}"
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
cd "$ROOT"
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
RS=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
TS=$(date -u +%Y%m%dT%H%M%SZ)
OUT="$H30/results/gate-$TS"
mkdir "$OUT"
S=scripts/analysis/hypotheses_2026q4
$RS $S/H30_foxo1_hotspots/01_export_roster_stage.R \
  RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/meta_matched.rds \
  "$OUT/roster_stage.tsv"
$PY -I $S/reads_extraction/h30_power_gate.py Analysis/Hypotheses/READS-20261009T141017Z/libraries "$OUT/gate_blind.tsv" 0.001 0.8
$PY -I $S/H30_foxo1_hotspots/02_gate_summary.py "$OUT/gate_blind.tsv" "$OUT/roster_stage.tsv" "$OUT/summary"
$PY -I -m pip freeze > "$OUT/pip_freeze.txt" 2>/dev/null || true
$RS -e 'writeLines(capture.output(sessionInfo()), file.path(commandArgs(TRUE)[1], "sessionInfo.txt"))' "$OUT"
echo "gate done: $OUT"
