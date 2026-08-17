#!/bin/bash
#SBATCH --job-name=ncvalidate
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/ncvalidate_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/ncvalidate_%j.err

set -euo pipefail
umask 027

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
PYTHON="/gpfs/commons/home/jameslee/micromamba/envs/alphagenome/bin/python"
SYSTEM_PYTHON="/usr/bin/python3"
RSCRIPT="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"
RUFF="/gpfs/commons/home/jameslee/.local/bin/ruff"

cd "${PROJECT_ROOT}"
"${SYSTEM_PYTHON}" -m pytest -q \
  scripts/manuscript/noncoding_resource/tests/test_audit_strandedness.py \
  scripts/manuscript/noncoding_resource/tests/test_gene_identity.py \
  scripts/manuscript/noncoding_resource/tests/test_annotate_program_lncrna_content.py \
  scripts/manuscript/noncoding_resource/tests/test_noncoding_genetics.py
"${PYTHON}" -m unittest \
  scripts/manuscript/noncoding_resource/tests/test_validate_program_without_lncrna.py -v
PYTHONPATH="${PROJECT_ROOT}/scripts/portal" \
  "${PYTHON}" -m unittest \
  scripts/portal/test_gene_catalog_v2.py \
  scripts/portal/test_gene_catalog_observability.py -v
"${RSCRIPT}" --vanilla scripts/manuscript/noncoding_resource/tests/test_bulk_lncrna_common.R
"${RSCRIPT}" --vanilla scripts/manuscript/noncoding_resource/tests/test_lncrna_biotype_association.R
"${RSCRIPT}" --vanilla -e \
  'parse(file="scripts/manuscript/noncoding_resource/build_stage_lncrna_candidate.R")'
"${RSCRIPT}" --vanilla -e \
  'parse(file="scripts/figures/noncoding_resource_fig3_panels.R")'
"${RSCRIPT}" --vanilla -e \
  'parse(file="scripts/figures/noncoding_resource_lncrna_supplement.R")'
"${RSCRIPT}" --vanilla -e \
  'parse(file="scripts/figures/noncoding_resource_plot_upgrade.R")'
"${RSCRIPT}" --vanilla -e \
  'parse(file="scripts/figures/noncoding_dna_plot_upgrade.R")'
for retired_script in \
  RNA-seq/53_ncrna_landscape.R \
  RNA-seq/55_ncrna_conservation.R \
  RNA-seq/56_ncrna_epigenomic.R \
  RNA-seq/57_ncrna_atlas_integration.R; do
  retired_log="$(mktemp)"
  if "${RSCRIPT}" --vanilla "${retired_script}" >"${retired_log}" 2>&1; then
    echo "Historical ncRNA script ran without explicit authorization: ${retired_script}" >&2
    rm -f "${retired_log}"
    exit 1
  fi
  grep -q "ALLOW_HISTORICAL_NCRNA_PIPELINE" "${retired_log}"
  rm -f "${retired_log}"
done
"${RUFF}" check \
  scripts/manuscript/noncoding_resource \
  scripts/portal/gene_catalog_v2.py \
  scripts/portal/test_gene_catalog_v2.py \
  scripts/portal/emit_gene_catalog_v2_contract.py \
  scripts/portal/validate_gene_catalog_v2_contract.py \
  scripts/portal/gene_catalog_observability.py \
  scripts/portal/test_gene_catalog_observability.py \
  scripts/figures/validate_noncoding_resource_plot_upgrade.py \
  scripts/manuscript/noncoding_resource/annotate_credible_set_vep.py \
  scripts/manuscript/noncoding_resource/build_noncoding_dna_candidate.py \
  scripts/manuscript/noncoding_resource/validate_noncoding_dna_candidate.py
"${RUFF}" format --check \
  scripts/manuscript/noncoding_resource \
  scripts/portal/gene_catalog_v2.py \
  scripts/portal/test_gene_catalog_v2.py \
  scripts/portal/emit_gene_catalog_v2_contract.py \
  scripts/portal/validate_gene_catalog_v2_contract.py \
  scripts/portal/gene_catalog_observability.py \
  scripts/portal/test_gene_catalog_observability.py \
  scripts/figures/validate_noncoding_resource_plot_upgrade.py \
  scripts/manuscript/noncoding_resource/annotate_credible_set_vep.py \
  scripts/manuscript/noncoding_resource/build_noncoding_dna_candidate.py \
  scripts/manuscript/noncoding_resource/validate_noncoding_dna_candidate.py
"${PYTHON}" scripts/manuscript/validate_resource_scope.py
