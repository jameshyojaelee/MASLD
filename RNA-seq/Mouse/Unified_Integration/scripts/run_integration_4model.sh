#!/bin/bash
#SBATCH --job-name=limma
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse/Unified_Integration/logs/integration_4model_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse/Unified_Integration/logs/integration_4model_%j.err
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=48:00:00
# NOTE: cpu partition (not bigmem) — bigmem enforces a 500G per-job minimum
# (QOSMinMemory); the mouse pipeline needs only ~128G, which fits the cpu cap (210G).

# =============================================================================
# run_integration_4model.sh — regenerate the mouse integration to the CANONICAL
# 4-diet-model roster (MCD, HFD, CDAHFD, FPC); excludes LIDPAD/GAN/AMLN.
#
# Differs from run_integration_slurm.sh in two ways:
#   1. SKIPS M02b (the Western/NASH/GSE246088 metabolic-overload pooling apparatus
#      — those datasets are excluded from config/mouse_datasets.yaml).
#   2. Archives any pre-existing per_diet/*_de_results.csv BEFORE M02 so that
#      M03 (which auto-globs per_diet/*_de_results.csv) meta-analyzes only the
#      4 freshly-written clean models — not the stale NASH_diet/Western files.
#
# Prereqs already applied in source (2026-06-16):
#   - M00 line 224 bug fixed (config$datasets -> mouse_cfg$datasets) so the
#     config dataset-exclusion filter actually runs.
#   - M02c DIETS reverted to c("MCD","HFD","CDAHFD","FPC").
# M02 needs no change (it is config-driven and splits CDAHFD_FPC -> CDAHFD+FPC).
#
# See docs/manuscript/NUMBERS.md "Mouse roster reconciliation (2026-06-16)".
# =============================================================================

# NOTE: activate the env BEFORE `set -u`. The conda binutils activation script
# references an unbound $ADDR2LINE, which `set -u` (nounset) turns into a fatal
# error and kills the job before any R runs (caused the 17299338 failure).
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

set -euo pipefail

BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$BASE/RNA-seq/Mouse/Unified_Integration"

SCRIPTS="scripts"
PERDIET="results/per_diet"

echo "=== Mouse Integration Pipeline (CANONICAL 4-MODEL roster) ==="
echo "Start: $(date)"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-NA}"
echo ""

# --- Archive any stale per-diet DE results so M03's glob sees only the fresh 4 ---
if compgen -G "$PERDIET/*_de_results.csv" > /dev/null || [ -f "$PERDIET/de_summary.csv" ]; then
  ARCH="$PERDIET/archive_pre_4model_${SLURM_JOB_ID:-manual}"
  mkdir -p "$ARCH"
  echo "Archiving stale per-diet results -> $ARCH"
  # only top-level CSVs; leave existing archive_*/ and backup subdirs untouched
  find "$PERDIET" -maxdepth 1 -type f -name "*.csv" -exec mv -v {} "$ARCH/" \;
  echo ""
fi

# M00: Harmonize metadata (now honors the 4-model config; bug fixed)
echo "=== M00: Harmonize Mouse Metadata ==="
Rscript "$SCRIPTS/M00_harmonize_mouse_metadata.R"
echo ""

# M01: Sample QC + count merge
echo "=== M01: Sample QC ==="
Rscript "$SCRIPTS/M01_mouse_sample_qc.R"
echo ""

# M02: Per-diet limma-voom DE (config-driven -> MCD, HFD, CDAHFD, FPC)
echo "=== M02: Per-Diet DE ==="
Rscript "$SCRIPTS/M02_mouse_per_diet_de.R"
echo ""

# M02b: SKIPPED (Western/NASH/GSE246088 pooling — excluded datasets).

# M02c: ashr shrinkage on the 4 per-diet results (DIETS reverted to 4-model)
echo "=== M02c: ashr shrinkage ==="
Rscript "$SCRIPTS/M02c_ashr_shrinkage.R"
echo ""

# M03: Meta-analysis (metafor rma across 4 diets) + pooled dream
echo "=== M03: Meta-Analysis + Pooled limma-voom-qw ==="
Rscript "$SCRIPTS/M03_mouse_meta_analysis.R"
echo ""

# M04: Consensus DEGs (atlas S2 tier input)
echo "=== M04: Consensus DEGs ==="
Rscript "$SCRIPTS/M04_mouse_consensus.R"
echo ""

echo "=== Mouse 4-model pipeline complete: $(date) ==="
echo "Sanity-check the new roster:"
echo "  awk -F, 'NR>1{print \$3}' metadata/unified_mouse_metadata.csv | sort | uniq -c   # diet_model"
echo "  ls $PERDIET/*_de_results.csv                                                      # expect MCD/HFD/CDAHFD/FPC only"
echo "  head -1 results/meta_analysis/meta_per_diet.csv"
