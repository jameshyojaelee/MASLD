#!/bin/bash
#SBATCH --job-name=H38_iga
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=logs/H38/%x-%j.out
#SBATCH --error=logs/H38/%x-%j.err
# H38: IgA fraction of heavy-chain reads. Usage: sbatch run_h38.sbatch <new_output_dir>  (H38_SMOKE=1 for the reduced run)
set -euo pipefail
OUT="${1:?new output dir required}"
BIN=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin
export PATH="$BIN:$PATH"
"$BIN/Rscript" scripts/analysis/hypotheses_2026q4/H38_iga/h38_iga.R "$OUT"
sha256sum scripts/analysis/hypotheses_2026q4/H38_iga/h38_iga.R > "$OUT/script_sha256.txt"
