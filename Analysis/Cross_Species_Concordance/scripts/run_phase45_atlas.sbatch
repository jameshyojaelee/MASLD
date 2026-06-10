#!/bin/bash
#SBATCH --job-name=xspecies
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Cross_Species_Concordance/logs/xspecies_p45_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Cross_Species_Concordance/logs/xspecies_p45_%j.err

# Re-run only Phases 4 (visualizations) + 5 (unified atlas) off the already-fresh
# Phase 1-3 C2 outputs. Phase 5 (the concordance_atlas_unified.csv deliverable) is run
# UNCONDITIONALLY even if Phase 4 (figures only) hiccups, so the atlas always regenerates.
set -uo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
set +u; micromamba activate rnaseq; set -u

SD="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Cross_Species_Concordance/scripts"
RES="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Cross_Species_Concordance/results"

echo "=== Phase 4 + 5 rerun (C2 atlas) ==="
echo "Date: $(date)  Node: $(hostname)"

echo "[04] Visualizations..."
Rscript "$SD/04_concordance_visualizations.R" || echo "WARN: Phase 4 (figures) returned non-zero; continuing to Phase 5."

echo "[05] Unified Concordance Atlas..."
Rscript "$SD/05_unified_concordance_atlas.R"
RC=$?

echo "=== DONE (Phase 5 rc=$RC) ==="
ls -la "$RES/concordance_atlas_unified.csv"
exit $RC
