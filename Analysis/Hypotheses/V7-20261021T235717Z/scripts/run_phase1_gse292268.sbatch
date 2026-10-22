#!/bin/bash
#SBATCH --job-name=V7p1_gse292268
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/V7-20261021T235717Z
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/V7-20261021T235717Z/logs/%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/V7-20261021T235717Z/logs/%j.err
# V7 phase 1 (gse292268): observability rows and in situ lineage of member genes.
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 PYTHONHASHSEED=0 OMP_NUM_THREADS=8
RUN=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/V7-20261021T235717Z
PY=/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python
sha256sum "$RUN/PRESPEC_V7_phase1.md" "$RUN/scripts/01_phase1.py" > "$RUN/script_sha256_phase1_gse292268.txt"
"$PY" -I "$RUN/scripts/01_phase1.py" gse292268
"$PY" -I -m pip freeze > "$RUN/pip_freeze_spatial.txt"
