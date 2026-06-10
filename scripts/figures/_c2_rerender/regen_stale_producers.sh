#!/bin/bash
#SBATCH --job-name=harmony
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/_c2_rerender/logs/regen_stale_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/_c2_rerender/logs/regen_stale_%j.err
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00

# ---------------------------------------------------------------------------
# Regenerate two stale/missing producer outputs that block figures, then
# re-render the dependent figures. All work runs on a compute node (cpu).
#
# OUTPUT 1: umap_coordinates.csv  (producer: 09_batch_correction_umap.R)
#   - STALE: predates the 5-cohort / 846-sample canonical merged_dge.rds.
#   - Re-run recomputes Harmony + UMAP on the current merged_dge.rds.
#   Dependent figs: fig1_umap.R, figS_batch_correction.R,
#                   figS_batch_correction_extra.R
#
# OUTPUT 2: conserved_enrichment.csv  (producer: 33_audit_sensitivity_analyses.R)
#   - MISSING at RNA-seq/results/audit_sensitivity/.
#   - Re-run reads concordance atlas + atlas/GSEA/DGIdb and writes Fisher
#     enrichment summary for Conserved_Core genes.
#   Dependent fig: figS_conserved_core.R (also reads umap_coordinates.csv)
# ---------------------------------------------------------------------------

# Activate env BEFORE `set -u` — the conda activate.d hooks reference unset vars
# (e.g. ADDR2LINE) and trip `set -u` if it is enabled first.
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

set -uo pipefail

BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export MASLD_PROJECT_ROOT="$BASE"
INT_SCRIPTS="$BASE/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts"
FIG="$BASE/scripts/figures"
UMAP_CSV="$BASE/RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/umap_coordinates.csv"
ENR_CSV="$BASE/RNA-seq/results/audit_sensitivity/conserved_enrichment.csv"
cd "$BASE"

echo "================================================================"
echo "[$(date '+%F %T')] START regen_stale_producers on $(hostname)"
echo "================================================================"

# ---------------------------------------------------------------------------
# PRODUCER 1: UMAP coordinates
# ---------------------------------------------------------------------------
echo ""
echo "[$(date '+%T')] >>> PRODUCER 1: 09_batch_correction_umap.R"
Rscript "$INT_SCRIPTS/09_batch_correction_umap.R"
rc1=$?
echo "[$(date '+%T')] <<< 09_batch_correction_umap.R rc=$rc1"
if [ -f "$UMAP_CSV" ]; then
  echo "UMAP CSV rows (incl header): $(wc -l < "$UMAP_CSV")"
  echo "UMAP CSV mtime: $(stat -c '%y' "$UMAP_CSV")"
fi

# ---------------------------------------------------------------------------
# PRODUCER 2: Conserved-core enrichment audit
# ---------------------------------------------------------------------------
echo ""
echo "[$(date '+%T')] >>> PRODUCER 2: 33_audit_sensitivity_analyses.R"
Rscript "$BASE/RNA-seq/33_audit_sensitivity_analyses.R"
rc2=$?
echo "[$(date '+%T')] <<< 33_audit_sensitivity_analyses.R rc=$rc2"
if [ -f "$ENR_CSV" ]; then
  echo "ENR CSV mtime: $(stat -c '%y' "$ENR_CSV")"
  echo "--- conserved_enrichment.csv ---"
  cat "$ENR_CSV"
  echo "--------------------------------"
fi

# ---------------------------------------------------------------------------
# DEPENDENT FIGURES
# ---------------------------------------------------------------------------
echo ""
echo "[$(date '+%T')] >>> RE-RENDER dependent figures"
RES="$FIG/_c2_rerender/results_stale_regen.tsv"
: > "$RES"
for name in fig1_umap figS_batch_correction figS_batch_correction_extra figS_conserved_core; do
  scr="$FIG/$name.R"
  log="$FIG/_c2_rerender/logs/${name}_staleregen.log"
  if [ ! -f "$scr" ]; then
    printf '%s\tMISSING\tscript file not found\n' "$name" >> "$RES"
    continue
  fi
  echo "[$(date '+%T')] render >>> $name"
  timeout 1200 Rscript "$scr" > "$log" 2>&1
  rc=$?
  if [ $rc -eq 0 ]; then
    printf '%s\tPASS\t\n' "$name" >> "$RES"
    echo "[$(date '+%T')] render <<< $name PASS"
  elif [ $rc -eq 124 ]; then
    printf '%s\tFAIL\tTIMEOUT after 1200s\n' "$name" >> "$RES"
    echo "[$(date '+%T')] render <<< $name TIMEOUT"
  else
    ferr=$(grep -m1 -iE "error|cannot|could not|no such file|not found|unable|object .* not found|subscript out of bounds" "$log" | head -1 | tr '\t' ' ' | cut -c1-300)
    [ -z "$ferr" ] && ferr=$(tail -1 "$log" | tr '\t' ' ' | cut -c1-300)
    printf '%s\tFAIL\t%s\n' "$name" "$ferr" >> "$RES"
    echo "[$(date '+%T')] render <<< $name FAIL rc=$rc :: $ferr"
  fi
done

echo ""
echo "=== RESULTS ($RES) ==="
cat "$RES"
echo ""
echo "[$(date '+%F %T')] DONE  producer1_rc=$rc1 producer2_rc=$rc2"
