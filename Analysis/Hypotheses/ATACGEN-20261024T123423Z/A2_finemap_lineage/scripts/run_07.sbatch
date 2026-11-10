#!/bin/bash
#SBATCH --job-name=A2_q2
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH -c 2
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/ATACGEN-20261024T123423Z/A2_finemap_lineage/logs/q2_%j.out
# Rerun of 07_q2 (adds the non-hepatocyte-only lineage breakdown) and 09_q4 into results/.
set -euo pipefail
A=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/ATACGEN-20261024T123423Z/A2_finemap_lineage
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 TMPDIR=$A/work
E=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin
for s in 07_q2 08_q3 09_q4; do echo "== $s $(date -u +%H:%M:%S)"; $E/python $A/scripts/$s.py; done
echo "== done $(date -u +%H:%M:%S)"
