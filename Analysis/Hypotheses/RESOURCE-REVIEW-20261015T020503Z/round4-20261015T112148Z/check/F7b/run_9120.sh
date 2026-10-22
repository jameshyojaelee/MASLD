#!/bin/bash
#SBATCH --job-name=9120
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/check/F7b/logs/%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/check/F7b/logs/%j.err
set -euo pipefail
D=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/check/F7b
export TMPDIR=$D/tmp XDG_CACHE_HOME=$D/tmp MPLCONFIGDIR=$D/tmp PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 OMP_NUM_THREADS=2
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
"$PY" -m pip freeze > $D/pip_freeze.txt 2>/dev/null || true
"$PY" $D/chk_r0802_r0803.py > $D/logs/chk_r0802_r0803.out
"$PY" $D/chk_r0808_r0809.py > $D/logs/chk_r0808_r0809.out
"$PY" $D/chk_r1101.py > $D/logs/chk_r1101.out
echo DONE
