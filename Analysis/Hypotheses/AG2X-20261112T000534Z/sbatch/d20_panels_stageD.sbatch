#!/bin/bash
#SBATCH --job-name=ag2x_d20_panels_stageD
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --time=90:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/AG2X-20261112T000534Z/logs/%x_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/AG2X-20261112T000534Z/logs/%x_%j.err
set -euo pipefail
O=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/AG2X-20261112T000534Z
PY=/gpfs/commons/home/jameslee/micromamba/envs/alphagenome2/bin/python
export PYTHONHASHSEED=20261111
cd $O/code
$PY -m pip freeze > $O/stageD/pip_freeze_d20_panels_stageD.txt
$PY -u d20_panels_stageD.py
$PY -u d21_captions_stageD.py
$PY -u d22_checks_stageD.py
for f in $O/stageD/panels/*/*.pdf; do ~/.local/bin/pdftoppm -png -r 300 -singlefile "$f" "${f%.pdf}"; done
