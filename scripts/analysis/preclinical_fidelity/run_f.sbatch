#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=90:00:00
# Model F: program scores (f1, R) then fidelity statistics (f2, Python). f2 runs once.
# Usage: sbatch --job-name=<n> --output=<exec>/logs/%j.out run_f.sbatch <exec_dir>
set -euo pipefail
repo=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mm=/gpfs/commons/home/jameslee/.local/bin/micromamba
pf=$repo/scripts/analysis/preclinical_fidelity
export PYTHONNOUSERSITE=1
[ -s "$1/scores/mouse_program_scores.tsv" ] || /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript "$pf/f1_score_programs.R" "$1/scores"
"$mm" run -n he_molecular python -m pip freeze > "$1/pip_freeze_f2.txt"
cd "$pf"
"$mm" run -n he_molecular python f2_fidelity.py --scores "$1/scores" --out "$1/f2_fidelity"
