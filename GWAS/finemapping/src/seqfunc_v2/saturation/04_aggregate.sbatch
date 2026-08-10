#!/usr/bin/env bash
#SBATCH --job-name=saturation
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=192G
#SBATCH --time=48:00:00
# Relative, resolved against the submission working directory, which
# submit_saturation_v2.sh sets to $SATURATION_V2_ROOT.  These were absolute
# paths into haplotype_saturation/v2/logs/, so a v3 saturation run would have
# written its logs into the v2 tree that the frozen release still references.
#SBATCH --output=logs/aggregate_%j.out
#SBATCH --error=logs/aggregate_%j.err
set -euo pipefail
ROOT=${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}
OUT=${SATURATION_V2_ROOT:-$ROOT/GWAS/finemapping/results/seqfunc/haplotype_saturation/v2}
mkdir -p "$OUT/logs"
python3 "$ROOT/GWAS/finemapping/src/seqfunc_v2/saturation/04_aggregate_scores.py"
