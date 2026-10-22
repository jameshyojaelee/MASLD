#!/bin/bash
#SBATCH --job-name=H23_verify
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H23-20261009T140416Z/logs/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H23-20261009T140416Z/logs/%x-%j.err
# Usage: sbatch run_h23_verify.sbatch <analysis_dir> <fits_dir> <new_output_dir> [n_refit] [comma-separated uids]
set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -u
Rscript scripts/analysis/hypotheses_2026q4/H23_multiancestry_finemapping/h23_verify_units.R "$1" "$2" "$3" "${4:-3}" ${5:+"$5"}
