#!/bin/bash
#SBATCH --job-name=N3_pseudobulk
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/N3-20261009T140326Z/logs/%x-%j.out
# N3 steps 1-2: target list from the saved Open Targets result, then donor-level hepatocyte pseudobulk.
# Usage: sbatch run_n3_pseudobulk.sbatch
set -euo pipefail
REPO=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
RES=$REPO/Analysis/Hypotheses/N3-20261009T140326Z
SCR=$REPO/scripts/analysis/hypotheses_2026q4/N3_target_expression
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
export PYTHONNOUSERSITE=1 PYTHONHASHSEED=0 OMP_NUM_THREADS=4
test -f "$RES/PRESPEC_N3.md"
(cd "$SCR" && sha256sum 01_build_target_list.py 02_hepatocyte_pseudobulk.py 03_donor_model.R) > "$RES/script_sha256.txt"
sha256sum "$RES/PRESPEC_N3.md" >> "$RES/script_sha256.txt"
"$PY" -m pip freeze > "$RES/pip_freeze.txt"
"$PY" -I "$SCR/01_build_target_list.py" "$RES"
"$PY" -I "$SCR/02_hepatocyte_pseudobulk.py" "$REPO" "$RES"
