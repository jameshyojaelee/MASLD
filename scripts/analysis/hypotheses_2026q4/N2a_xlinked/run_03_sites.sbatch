#!/bin/bash
#SBATCH --job-name=N2a_sites
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/N2a-20261009T202200Z/logs/slurm-%x-%j.out
# Builds inputs/sites.tsv (common exonic SNVs on chrX non-PAR and chr19-22). Reads no BAM.
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/N2a-20261009T202200Z
S=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/hypotheses_2026q4/N2a_xlinked
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
sha256sum "$S/03_build_site_lists.py" > "$ROOT/logs/script_sha256_03.txt"
"$PY" -I "$S/03_build_site_lists.py" "$ROOT"
