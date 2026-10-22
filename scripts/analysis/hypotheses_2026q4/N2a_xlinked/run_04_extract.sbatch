#!/bin/bash
#SBATCH --job-name=N2a_extract
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/N2a-20261009T202200Z/logs/slurm-%x-%j.out
# Allele counts at the N2a sites for the 64 libraries of the included pairs. Single job, no array.
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/N2a-20261009T202200Z
S=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/hypotheses_2026q4/N2a_xlinked
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
sha256sum "$S/04_extract_counts.py" "$S/n2a_lib.py" >> "$ROOT/logs/script_sha256_04.txt"
"$PY" -m pip freeze > "$ROOT/logs/pip_freeze_extract.txt"
"$PY" -I "$S/04_extract_counts.py" "$ROOT" 15
