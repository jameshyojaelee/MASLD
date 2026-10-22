#!/bin/bash
#SBATCH --job-name=V5bridge
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/V5-20261021T122723Z
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/V5-20261021T122723Z/logs/%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/V5-20261021T122723Z/logs/%j.err
# V5 job 2: discovery bridge (G4).
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 PYTHONHASHSEED=0 OMP_NUM_THREADS=4
RUN=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/V5-20261021T122723Z
H90S=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H90-20261015T120337Z/scripts
sha256sum "$RUN/PRESPEC_V5.md" "$RUN"/scripts/*.py "$H90S/h90_common.py" > "$RUN/script_sha256_job2.txt"
/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python -I "$RUN/scripts/02_bridge.py"
