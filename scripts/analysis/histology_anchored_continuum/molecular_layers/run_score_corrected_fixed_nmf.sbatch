#!/bin/bash
#SBATCH --job-name=fixed_nmf_score
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/logs/fixed_nmf_corrected/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/logs/fixed_nmf_corrected/%x-%j.err

# Corrected-count fixed-W NMF scores for the 774 primary-arm units of the
# bulk-corrected-controls-20261008T162308Z refit, under the r2 transport protocol.
# Bases unchanged; unavailable axes (k4-P2, k6-P6 for the stage variant) carry NA.
set -euo pipefail
umask 027
project_root="/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
env_bin="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin"
cand="${project_root}/figures/candidates/resource-corrections-20261005-v1"
fit="${project_root}/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z"
ts="$(date -u +%Y%m%dT%H%M%SZ)"
out="${project_root}/RNA-seq/results/manuscript_release/candidates/fixed-nmf-corrected-scores-${ts}"
test ! -e "${out}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"${env_bin}/Rscript" --vanilla \
  "${project_root}/scripts/analysis/histology_anchored_continuum/molecular_layers/score_corrected_fixed_nmf.R" \
  "${cand}/bulk_execution_sources/revised_fixed_nmf_projection_r7_review/revised_fixed_nmf_projection_protocol_with_transport_checks_r2.rds" \
  "${cand}/verified_counts_r8_20261007/assembled_774_primary" \
  "${cand}/identity_resource_qualified/bulk_execution_inputs/primary_bulk_metadata.tsv" \
  "${fit}/primary_source_controls/model_input_manifest.tsv" \
  "${out}"
echo "FIXED_NMF_CORRECTED_SCORES_COMPLETE ${out}"
