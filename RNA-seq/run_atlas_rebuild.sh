#!/bin/bash
#SBATCH --job-name=atlas_rebuild
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/atlas_rebuild_%j.log
#SBATCH --error=logs/atlas_rebuild_%j.log

# Full atlas rebuild (T0.2 update 2026-04-22):
#   27a -> L8 ATAC -> 40b (sc-TWAS) -> 45a (spatial + scRNA + LIANA + sources_active)
#   -> 75 (causal overhaul: HyPrColoc / cTWAS / SuSiE merge)
#   -> 57 (ncRNA integration: ceRNA + synteny + Elatus sc-lncRNA + classification)
#   -> 217 (stratified causal: sex/subtype/progression/spatial/pharma merges)
#   -> 27b (benchmark presets)
# CLAUDE.md canonical execution order: 27a -> 75 -> 217 (+57 before 217 for ncRNA cols).
# If you need the stratified-causal axes to be non-stale you must also rerun
# 207/210/211/213/215 before 217. Those are heavy and are submitted separately via
# docs/archive/code_review_2026-04-21/fix_progress/run_atlas_*.sh (see T0.2 fix chain;
# the code_review tree was archived on 2026-05-24, the .sh scripts remain executable at the archived path).

set -eo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
mkdir -p logs

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

echo "=== Atlas Rebuild ==="
echo "Job ID: ${SLURM_JOB_ID:-interactive}"
echo "Start: $(date)"

echo ""
echo "--- Step 1: 27a (core assembly + L8 ATAC + ferroptosis/zonation) ---"
Rscript 27a_assemble_evidence_atlas.R

echo ""
echo "--- Step 1b: 35_atac_integration (L8 ATAC/SCENIC+/chromVAR) ---"
echo "NOTE: Requires snapatac2 env. HARD step (handoff 2026-05-31): NO swallow — a Script-35 failure aborts"
echo "      the rebuild rather than silently producing a 0-column L8 layer."
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
(micromamba activate snapatac2 && python Analysis/ATAC/Integration/scripts/35_atac_integration.py)
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq

echo ""
echo "--- Step 2: 40b (sc-TWAS integration into L4) ---"
# 40b_integrate_sceqtl_twas.R was archived in the MR ditch (2026-04-22) — skip if absent so the
# rebuild does not fatally abort under `set -eo pipefail` (the L8 ATAC step uses the same pattern).
if [ -f 40b_integrate_sceqtl_twas.R ]; then
  Rscript 40b_integrate_sceqtl_twas.R
else
  echo "WARNING: 40b_integrate_sceqtl_twas.R not found (archived 2026-04-22) — skipping sc-TWAS step."
fi

echo ""
echo "--- Step 3: 45a (spatial + scRNA + LIANA + sources_active recompute) ---"
Rscript 45a_integrate_all_sources.R

echo ""
echo "--- Step 4: 75 (causal overhaul: SuSiE + HyPrColoc + cTWAS) ---"
Rscript 75_integrate_causal_overhaul.R

echo ""
echo "--- Step 5: 57 (ncRNA atlas integration: ceRNA + synteny + sc-lncRNA) ---"
Rscript 57_ncrna_atlas_integration.R

echo ""
echo "--- Step 6: 217 (stratified causal atlas integration) ---"
Rscript 217_stratified_causal_atlas.R

echo ""
echo "--- Step 7: 27b (benchmark presets) ---"
Rscript 27b_benchmark_presets.R

echo ""
echo "--- Assert: L8 ATAC layer landed in the atlas (handoff 2026-05-31; was silently 0) ---"
Rscript -e 'h <- names(data.table::fread("results/multi_evidence/multi_evidence_atlas.csv", nrows=0));
  l8 <- grep("hepatocyte_da_|scenic_regulon|chromvar|^atac_|disease_regulon", h, value=TRUE);
  cat("  L8 columns present:", length(l8), "->", paste(head(l8,15), collapse=", "), "\n");
  if (length(l8) < 1) stop("L8 ATAC layer ABSENT from atlas — Script 35 / 45a integration failed (do NOT ship a 0-L8 atlas)")'

echo ""
echo "=== Atlas rebuild complete ==="
echo "Check: RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"
echo "End: $(date)"
