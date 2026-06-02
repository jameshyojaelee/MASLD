#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=4:00:00
#SBATCH --job-name=phase10_sc
#SBATCH --output=Analysis/SingleCell/logs/phase10_sc_%j.out
#SBATCH --error=Analysis/SingleCell/logs/phase10_sc_%j.err

set -eo pipefail

# Single-cell scripts here are R-based and depend only on canonical
# gene_level_coloc.csv (not the atlas), so they run in parallel with the
# atlas chain. Use rnaseq env (R + data.table; no GPU dependencies).
source /gpfs/commons/home/jameslee/.bashrc
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

mkdir -p Analysis/SingleCell/logs

PHASE_LOG="Analysis/SingleCell/logs/phase10_sc_phase_timings_${SLURM_JOB_ID:-$$}.log"
echo "Phase 10 single-cell COLOC enrichment — started $(date)" | tee "$PHASE_LOG"

run_phase () {
  local label="$1"; shift
  local cmd="$*"
  local t0=$(date +%s)
  echo "[$(date)] BEGIN: $label" | tee -a "$PHASE_LOG"
  echo "  CMD: $cmd" | tee -a "$PHASE_LOG"
  if eval "$cmd"; then
    local t1=$(date +%s)
    echo "[$(date)] PASS: $label  ($((t1 - t0))s)" | tee -a "$PHASE_LOG"
  else
    local rc=$?
    local t1=$(date +%s)
    echo "[$(date)] FAIL: $label  ($((t1 - t0))s)  rc=$rc" | tee -a "$PHASE_LOG"
    return $rc
  fi
}

run_phase "01_312c_hepatocyte_enrichment"  "Rscript Analysis/SingleCell/scripts/312c_hepatocyte_enrichment.R"  || echo "WARN: 312c failed, continuing"
run_phase "02_483_mcp_coloc_enrichment"    "Rscript Analysis/SingleCell/scripts/400_mcp/483_coloc_enrichment.R" || echo "WARN: 483 failed, continuing"

echo "[$(date)] PHASE 10 SINGLE-CELL COMPLETE" | tee -a "$PHASE_LOG"
