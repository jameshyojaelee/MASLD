#!/bin/bash
#SBATCH --job-name=p2a_atlas
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/multi_evidence/logs/p2a_atlas_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/multi_evidence/logs/p2a_atlas_%j.err

# Phase 2A: rebuild multi-evidence atlas under 5-cohort dream.
#
# Chain order (per CLAUDE.md + existing run_atlas_rebuild.sh):
#   27a  core assembly + L8 ATAC merge + ferroptosis/zonation
#   35   L8 ATAC integration  (snapatac2 env; skip if unavailable)
#   40b  sc-TWAS integration
#   45a  spatial + scRNA + LIANA + sources_active recompute
#   75   causal overhaul (HyPrColoc + cTWAS + SuSiE merge)
#   57   ncRNA atlas integration (ceRNA + synteny + sc-lncRNA)
#   217  stratified causal atlas merges
#   27b  benchmark presets

set -eo pipefail

BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p "$BASE/RNA-seq/results/multi_evidence/logs"
cd "$BASE/RNA-seq"

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

echo "============================================================"
echo "  Phase 2A — atlas rebuild (5-cohort)"
echo "  Started: $(date)"
echo "============================================================"

echo
echo "--- 27a: core assembly ---"
Rscript 27a_assemble_evidence_atlas.R

echo
echo "--- 35: ATAC integration (snapatac2 env, may skip) ---"
cd "$BASE"
if (micromamba activate snapatac2 && python Analysis/ATAC/Integration/scripts/35_atac_integration.py); then
  echo "  L8 ATAC integration OK"
else
  echo "  WARNING: L8 ATAC integration skipped (env or data not available). S5 epigenomic columns may be 0%."
fi
micromamba activate rnaseq
cd "$BASE/RNA-seq"

echo
echo "--- 40b: sc-TWAS integration ---"
Rscript 40b_integrate_sceqtl_twas.R

echo
echo "--- 45a: spatial + scRNA + LIANA + sources_active ---"
Rscript 45a_integrate_all_sources.R

echo
echo "--- 75: causal overhaul ---"
Rscript 75_integrate_causal_overhaul.R

echo
echo "--- 57: ncRNA atlas integration ---"
Rscript 57_ncrna_atlas_integration.R

echo
echo "--- 217: stratified causal atlas integration ---"
Rscript 217_stratified_causal_atlas.R

echo
echo "--- 27b: benchmark presets ---"
Rscript 27b_benchmark_presets.R

echo
echo "============================================================"
echo "  Phase 2A complete: $(date)"
echo "  Atlas: $BASE/RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"
echo "============================================================"
