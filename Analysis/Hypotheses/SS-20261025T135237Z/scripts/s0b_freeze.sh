#!/bin/bash
# SS job S0b: copy the frozen discovery references used by the SteatoSITE replication into frozen/, with SHA256SUMS.
# Run before any SteatoSITE expression is read. Sources are only read.
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
SS=$ROOT/Analysis/Hypotheses/SS-20261025T135237Z
F=$SS/frozen
[ -n "$(ls -A $F 2>/dev/null)" ] && { echo "frozen/ not empty; refusing"; exit 1; }
D=$ROOT/Analysis/Hypotheses/DIG1-20261024T113951Z
H51=$ROOT/Analysis/Hypotheses/H51-20261015T112617Z
cp_one() { mkdir -p "$F/$2"; cp -p "$1" "$F/$2/"; echo -e "$2/$(basename $1)\t$1" >> "$F/SOURCES.tsv"; }
cp_one $H51/stage1/D0_discovery_grouped_per_gene.tsv.gz h51
cp_one $H51/stage2-20261015T164829Z/per_gene_states.tsv.gz h51
cp_one $H51/PRESPEC_H51.md h51
cp_one $H51/RESULT_H51.md h51
cp_one $ROOT/figures/candidates/bulk-release-corrected-controls-20261008T162308Z/analysis/stage_extensions/stage_extension_all_gene_results.tsv stage
cp_one $ROOT/figures/main/fig3_bulk_transcriptomics/panels/data/panel_deep_20261017/fig3e.tsv stage
cp_one $ROOT/figures/main/fig3_bulk_transcriptomics/panels/data/panel_deep_20261017/fig3f.tsv stage
cp_one $ROOT/figures/main/fig3_bulk_transcriptomics/manifests/fig3f_corrected_20261005.tsv stage
for f in transcript_stage_vs_continuum.tsv.gz transcript_stage_trend_cohort.tsv.gz hotspot_stage_vs_continuum.tsv; do
  cp_one $ROOT/figures/candidates/continuum-774-panels-20261011T165944Z/stage_continuum_comparison/source_tables/$f stage_continuum; done
cp_one $ROOT/figures/candidates/continuum-774-panels-20261011T165944Z/stage_continuum_comparison/comparison_summary.tsv stage_continuum
cp_one $ROOT/Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/program_registry_v2.tsv programs
cp_one $ROOT/RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv release
for f in out_genotype_state/primary_pnpla3.tsv out_genotype_state/secondary_holm.tsv out_panx1_ase/D_dosage_models.tsv \
         out_release_vs_disease_loci/loci.tsv out_release_vs_disease_loci/coloc_genes_lead_z.tsv \
         out_genome_disease_state/gene_disease_state.tsv.gz out_genotype_table/site_qc.tsv scripts/62_targets.tsv \
         scripts/62_targets_annotated.tsv scripts/71_targets.tsv; do cp_one $D/$f dig1; done
cp_one $ROOT/figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/passport_gene_index.parquet catalog
cd $F && find . -type f ! -name SHA256SUMS ! -name SOURCES.tsv | sort | xargs sha256sum > SHA256SUMS
wc -l SHA256SUMS
