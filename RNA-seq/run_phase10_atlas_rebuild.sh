#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=8:00:00
#SBATCH --job-name=phase10_atlas
#SBATCH --output=RNA-seq/logs/phase10_atlas_%j.out
#SBATCH --error=RNA-seq/logs/phase10_atlas_%j.err

set -eo pipefail

# Activate rnaseq env
source /gpfs/commons/home/jameslee/.bashrc
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

mkdir -p RNA-seq/logs

PHASE_LOG="RNA-seq/logs/phase10_atlas_phase_timings_${SLURM_JOB_ID:-$$}.log"
echo "Phase 10 atlas rebuild — started $(date)" | tee "$PHASE_LOG"

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

# =====================================================================
# Step 1: Atlas reassembly + downstream causal integration
# =====================================================================
run_phase "01_27a_assemble_atlas"          "Rscript RNA-seq/27a_assemble_evidence_atlas.R"
run_phase "02_75_integrate_causal"         "Rscript RNA-seq/75_integrate_causal_overhaul.R"
run_phase "03_48_cross_ancestry"           "Rscript RNA-seq/48_cross_ancestry_replication.R"
run_phase "04_202_geneset"                 "Rscript RNA-seq/202_geneset_enrichment.R"
run_phase "05_205_proportion_coloc"        "Rscript RNA-seq/205_proportion_coloc.R"
run_phase "06_200_intact"                  "Rscript RNA-seq/200_intact_scoring.R"
run_phase "07_58_diamond_kda"              "Rscript RNA-seq/58_diamond_kda.R"

# =====================================================================
# Step 2: Stratified causal projections (sequential)
# =====================================================================
run_phase "08_207_sex"                     "Rscript RNA-seq/207_sex_stratified_coloc.R"           || echo "WARN: 207 failed, continuing"
run_phase "09_208_subtype"                 "Rscript RNA-seq/208_subtype_coloc.R"                  || echo "WARN: 208 failed, continuing"
run_phase "10_210_progression"             "Rscript RNA-seq/210_progression_coloc.R"              || echo "WARN: 210 failed, continuing"
run_phase "11_211_progression_driver"      "Rscript RNA-seq/211_progression_driver_genetics.R"    || echo "WARN: 211 failed, continuing"
run_phase "12_213_spatial"                 "Rscript RNA-seq/213_spatial_coloc.R"                  || echo "WARN: 213 failed, continuing"
run_phase "13_215_pharmacogenomic"         "Rscript RNA-seq/215_pharmacogenomic_targets.R"        || echo "WARN: 215 failed, continuing"

# =====================================================================
# Step 3: Stratified causal merge + Bayesian convergence
# =====================================================================
run_phase "14_217_stratified_atlas"        "Rscript RNA-seq/217_stratified_causal_atlas.R"
run_phase "15_46d_convergence_evidence"    "Rscript RNA-seq/46d_convergence_evidence.R"

echo "[$(date)] PHASE 10 ATLAS REBUILD COMPLETE" | tee -a "$PHASE_LOG"
