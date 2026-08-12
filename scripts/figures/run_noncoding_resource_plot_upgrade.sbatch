#!/bin/bash
#SBATCH --job-name=ncplots
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/ncplots_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/ncplots_%j.err

set -euo pipefail
umask 027

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RNASEQ_ENV="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq"
NC_ROOT="${PROJECT_ROOT}/RNA-seq/results/noncoding_resource/candidates/noncoding-bulk-2026-08-11"
PROGRAM_ROOT="${PROJECT_ROOT}/Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot"
CATALOG_ROOT="${PROJECT_ROOT}/Analysis/Multimodal_Program_Projection/candidates/gene-catalog-v2-contract-2026-08-11"
CATALOG_REVIEW="${PROJECT_ROOT}/Analysis/Multimodal_Program_Projection/candidates/gene-catalog-v2-contract-review-2026-08-11/VALIDATED.json"
OUTPUT="${PROJECT_ROOT}/figures/candidates/noncoding-resource-2026-08-12/plot-upgrade-v3"

test ! -e "${OUTPUT}"
mkdir -p "$(dirname "${OUTPUT}")"
export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export PATH="${RNASEQ_ENV}/bin:${PATH}"

"${RNASEQ_ENV}/bin/Rscript" --vanilla scripts/figures/noncoding_resource_plot_upgrade.R \
  --bulk-counts="${NC_ROOT}/biotype_association/run_gate19751171/bulk_deg_counts_by_biotype.tsv" \
  --biotype-association="${NC_ROOT}/biotype_association/run_gate19751171/biotype_deg_association.tsv" \
  --biotype-review="${NC_ROOT}/biotype_association/review_gate19751171_job19751607/VALIDATED.json" \
  --stage-counts="${NC_ROOT}/stage/run_gate19751171/stage_counts_by_biotype.tsv" \
  --stage-review="${NC_ROOT}/stage/run_gate19751171.review/VALIDATED.json" \
  --bulk="${NC_ROOT}/bulk_disease/run_gate19751171/bulk_lncrna_results.tsv" \
  --cohort="${NC_ROOT}/bulk_disease/run_gate19751171/bulk_lncrna_cohort_effects.tsv" \
  --loco="${NC_ROOT}/bulk_disease/run_gate19751171/bulk_lncrna_leave_one_cohort_out.tsv" \
  --bulk-review="${NC_ROOT}/bulk_disease/review_gate19751171/VALIDATED.json" \
  --program-content="${NC_ROOT}/program_annotation/run_19749485/program_lncrna_content.tsv" \
  --program-sensitivity="${NC_ROOT}/program_without_lncrna/combined/program_without_lncrna_sensitivity.tsv" \
  --program-review="${NC_ROOT}/program_without_lncrna/combined/VALIDATED.json" \
  --program-registry="${PROGRAM_ROOT}/program_registry_v2.tsv" \
  --catalog-rules="${CATALOG_ROOT}/masld_gene_catalog_v2_experiment_rules.tsv" \
  --catalog-review="${CATALOG_REVIEW}" \
  --output="${OUTPUT}"
