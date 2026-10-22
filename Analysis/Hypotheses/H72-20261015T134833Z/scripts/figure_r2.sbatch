#!/bin/bash
#SBATCH --job-name=6125
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H72-20261015T134833Z/logs/6125-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H72-20261015T134833Z/logs/6125-%j.err
# Cosmetic re-render of the H72 panels (right margin, two-row legend) into a new leaf; same data, same marks.
set -euo pipefail
H72=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H72-20261015T134833Z
OUT=${H72}/run-20261015T142309Z
FIG=${OUT}/figures-r2-$(date -u +%Y%m%dT%H%M%SZ)
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript ${H72}/scripts/02_h72_figure.R "${OUT}" "${FIG}" 2> "${FIG}.captions.txt"
