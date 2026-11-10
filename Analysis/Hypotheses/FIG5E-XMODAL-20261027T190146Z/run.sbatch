#!/bin/bash
#SBATCH --job-name=fig5e-xmodal
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=Analysis/Hypotheses/FIG5E-XMODAL-20261027T190146Z/logs/%j.out
#SBATCH --error=Analysis/Hypotheses/FIG5E-XMODAL-20261027T190146Z/logs/%j.err
set -euo pipefail
P=/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python
$P Analysis/Hypotheses/FIG5E-XMODAL-20261027T190146Z/xmodal.py Analysis/Hypotheses/FIG5E-XMODAL-20261027T190146Z/out
$P -m pip freeze > Analysis/Hypotheses/FIG5E-XMODAL-20261027T190146Z/out/pip_freeze.txt
