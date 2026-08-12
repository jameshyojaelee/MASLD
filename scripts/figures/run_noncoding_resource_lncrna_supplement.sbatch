#!/bin/bash
#SBATCH --job-name=nclnfig
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/nclnfig_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/nclnfig_%j.err

set -euo pipefail
umask 027

ROOT="RNA-seq/results/noncoding_resource/candidates/noncoding-bulk-2026-08-11/bulk_disease"
OUT="figures/candidates/noncoding-resource-2026-08-11/lncrna_supplement_2026-08-11-r2"
export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"

/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript --vanilla \
  scripts/figures/noncoding_resource_lncrna_supplement.R \
  --bulk="${ROOT}/run_gate19751171/bulk_lncrna_results.tsv" \
  --cohort="${ROOT}/run_gate19751171/bulk_lncrna_cohort_effects.tsv" \
  --loco="${ROOT}/run_gate19751171/bulk_lncrna_leave_one_cohort_out.tsv" \
  --review="${ROOT}/review_gate19751171/VALIDATED.json" \
  --output="${OUT}"
