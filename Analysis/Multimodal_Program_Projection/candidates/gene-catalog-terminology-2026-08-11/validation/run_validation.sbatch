#!/usr/bin/env bash
#SBATCH --job-name=catalog
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=Analysis/Multimodal_Program_Projection/candidates/gene-catalog-terminology-2026-08-11/validation/slurm-%j.out
#SBATCH --error=Analysis/Multimodal_Program_Projection/candidates/gene-catalog-terminology-2026-08-11/validation/slurm-%j.err

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

python3 scripts/manuscript/validate_resource_scope.py
python3 -m py_compile \
  scripts/manuscript/validate_resource_scope.py \
  scripts/portal/build_evidence_passport_ui.py \
  scripts/portal/build_evidence_passport_bundle.py \
  scripts/portal/generate_evidence_passports.py \
  scripts/portal/validate_evidence_passports.py \
  scripts/portal/prepare_passport_attestations.py \
  scripts/portal/test_evidence_passport_bundle.py \
  scripts/portal/stage_hf_dataset.py \
  scripts/portal/_portal_io.py \
  scripts/figures/fig5_gene_catalog_panels.py \
  scripts/manuscript/program_context_v2/coordinator/build_scientific_validation_registry.py \
  scripts/manuscript/program_context_v2/coordinator/coordinator_contract.py \
  scripts/manuscript/program_context_v2/coordinator/prepare_real_release.py \
  scripts/manuscript/program_context_v2/rel04_build_candidate_manuscript.py \
  scripts/manuscript/program_context_v2/release_common.py \
  scripts/manuscript/program_context_v2/release_products.py \
  scripts/manuscript/program_context_v2/tests/fixture_factory.py \
  masld-atlas-v2/scripts/preprocess_atlas_data.py \
  streamlit_convergence/app.py
python3 -m json.tool config/resource_paper_scope.config.json >/dev/null

bash -n \
  scripts/figures/run_fig5_gene_catalog_panels.sbatch \
  scripts/portal/run_gene_catalog.sbatch \
  scripts/portal/run_gene_catalog_ui_smoke.sbatch \
  scripts/portal/run_gene_catalog_ui_tests.sbatch \
  scripts/portal/rebuild_web_data.sbatch \
  RNA-seq/run_atlas_rebuild_chain.sh

micromamba run -n rnaseq Rscript -e \
  'parse(file="scripts/figures/fig5_gene_catalog_panels.R"); parse(file="scripts/figures/load_figure_data.R")'
node --check masld-atlas-v2/scripts/check-resource-release-gate.mjs

cd masld-atlas-v2
npx eslint \
  src/app/layout.tsx \
  src/app/downloads/page.tsx \
  src/app/drugs/page.tsx \
  src/app/genetics/genetics-client.tsx \
  src/app/genetics/page.tsx \
  src/app/programs/page.tsx \
  src/app/proteomics/page.tsx \
  src/components/hero/hero.tsx \
  src/components/sidebar.tsx \
  src/lib/atlas-constants.ts

printf 'GENE_CATALOG_VALIDATION\tPASS\n'
