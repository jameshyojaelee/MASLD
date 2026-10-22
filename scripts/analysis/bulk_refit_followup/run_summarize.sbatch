#!/bin/bash
#SBATCH --job-name=refit_followup_summary
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/logs/bulk_refit_followup/%x_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/logs/bulk_refit_followup/%x_%j.err

# Summaries for the leave-one-cohort-out refits and the composition rebuild.
#   sbatch --job-name=<name> --export=ALL,SUMMARY_SCRIPT=summarize_loco.R,LOCO_TS=...,LOCO_SUMMARY_ROOT=... run_summarize.sbatch
#   sbatch --job-name=<name> --export=ALL,SUMMARY_SCRIPT=summarize_composition.R,COMP_ROOT=...,COMP_SUMMARY_ROOT=... run_summarize.sbatch
set -euo pipefail
umask 027
: "${SUMMARY_SCRIPT:?SUMMARY_SCRIPT is required}"
rnaseq_env="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq"
export LANG="en_US.UTF-8" LC_ALL="en_US.UTF-8"
export PATH="${rnaseq_env}/bin:${PATH}"
"${rnaseq_env}/bin/Rscript" --vanilla \
  "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/bulk_refit_followup/${SUMMARY_SCRIPT}"
