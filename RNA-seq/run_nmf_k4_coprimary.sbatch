#!/bin/bash
#SBATCH --partition=bigmem --qos=nslab
#SBATCH --cpus-per-task=16 --mem=500G --time=48:00:00
#SBATCH --job-name=nmf_k4_coprimary
#SBATCH --output=RNA-seq/logs/nmf_k4_coprimary_%j.out
#SBATCH --error=RNA-seq/logs/nmf_k4_coprimary_%j.err
# =============================================================================
# run_nmf_k4_coprimary.sbatch — NMF k=4 co-primary analysis
#
# Runs three steps sequentially:
#   Step 1: 95b seed sensitivity at k=4 (3 seeds x nrun=50)
#           -> nmf_k4_seed_robustness.csv
#   Step 2: Script 44 with NMF_CHOSEN_K=4 (full downstream: assignments,
#           markers, pathways, clinical, figures)
#           -> nmf_k4_programs.csv (renamed from nmf_assignments.csv)
#   Step 3: Consensus k-means ARI sanity check
#           -> nmf_vs_kmeans_ari.csv
#
# Uses existing clean cache (k=4 fit already present).
# =============================================================================
set -eo pipefail
mkdir -p RNA-seq/logs
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

OUT_DIR="RNA-seq/results/subtypes"

echo "============================================================"
echo "  NMF k=4 co-primary pipeline"
echo "  Start: $(date)"
echo "  SLURM_JOB_ID: ${SLURM_JOB_ID}"
echo "  SLURM_CPUS_PER_TASK: ${SLURM_CPUS_PER_TASK}"
echo "============================================================"

# ---- Step 1: Seed sensitivity (3 seeds x nrun=50 x k=4) ----
echo ""
echo "=== Step 1/3: Seed sensitivity at k=4 ==="

export NMF_K_LIST="4"
export NMF_RUNS=50
export NMF_REUSE_SEED42_CACHE=FALSE
export NMF_CACHE_PATH="RNA-seq/results/subtypes/nmf_results_cache_clean.rds"

Rscript RNA-seq/95b_nmf_seed_sensitivity.R

# The default output from 95b overwrites nmf_ksweep_seed_robustness.csv with
# only k=4 rows. We want to preserve the existing k=6,8 data and write k=4
# separately. Post-process: extract k=4 rows and restore original.
echo "--- Post-processing seed robustness output ---"
if [ -f "${OUT_DIR}/nmf_ksweep_seed_robustness.csv" ]; then
  # 95b writes k=4 rows to the main file. Since we want a separate k=4 file,
  # we need to handle this carefully.
  # Strategy: 95b overwrites the file. Extract k=4 rows to separate file,
  # then restore original (with k=6,8 data) + append k=4.

  # The current file now has ONLY k=4 rows (95b overwrites).
  # Save as k=4 specific file.
  cp "${OUT_DIR}/nmf_ksweep_seed_robustness.csv" \
     "${OUT_DIR}/nmf_k4_seed_robustness.csv"
  echo "Wrote: ${OUT_DIR}/nmf_k4_seed_robustness.csv"

  # Restore from git if the original had k=6,8 data, then append k=4
  if git show HEAD:"${OUT_DIR}/nmf_ksweep_seed_robustness.csv" > /dev/null 2>&1; then
    git show HEAD:"${OUT_DIR}/nmf_ksweep_seed_robustness.csv" > "${OUT_DIR}/nmf_ksweep_seed_robustness.csv.bak"
    # Combine: header + k=6,8 body + k=4 body (skip header)
    head -1 "${OUT_DIR}/nmf_ksweep_seed_robustness.csv.bak" > "${OUT_DIR}/nmf_ksweep_seed_robustness.csv"
    tail -n +2 "${OUT_DIR}/nmf_ksweep_seed_robustness.csv.bak" >> "${OUT_DIR}/nmf_ksweep_seed_robustness.csv"
    tail -n +2 "${OUT_DIR}/nmf_k4_seed_robustness.csv" >> "${OUT_DIR}/nmf_ksweep_seed_robustness.csv"
    rm -f "${OUT_DIR}/nmf_ksweep_seed_robustness.csv.bak"
    echo "Restored combined seed robustness (k=4 + k=6,8)"
  fi
fi

echo ""
echo "=== Step 1 complete ==="

# ---- Step 2: Full k=4 downstream via Script 44 ----
echo ""
echo "=== Step 2/3: Script 44 at k=4 ==="

# Back up existing k=6 outputs (if present) before Script 44 overwrites them.
for f in nmf_assignments.csv program_labels.csv subtype_markers.csv \
         subtype_pathways.csv subtype_clinical_associations.csv; do
  if [ -f "${OUT_DIR}/${f}" ]; then
    cp "${OUT_DIR}/${f}" "${OUT_DIR}/${f}.k6_backup"
  fi
done

export NMF_CHOSEN_K=4
export NMF_CACHE_PATH="RNA-seq/results/subtypes/nmf_results_cache_clean.rds"
unset NMF_K_LIST NMF_RUNS NMF_REUSE_SEED42_CACHE

Rscript RNA-seq/44_molecular_subtyping.R

# Rename k=4 outputs to k4-specific names
echo "--- Renaming k=4 outputs ---"
cp "${OUT_DIR}/nmf_assignments.csv" "${OUT_DIR}/nmf_k4_programs.csv"
cp "${OUT_DIR}/program_labels.csv" "${OUT_DIR}/nmf_k4_program_labels.csv"
cp "${OUT_DIR}/subtype_markers.csv" "${OUT_DIR}/nmf_k4_subtype_markers.csv"
cp "${OUT_DIR}/subtype_pathways.csv" "${OUT_DIR}/nmf_k4_subtype_pathways.csv"
cp "${OUT_DIR}/subtype_clinical_associations.csv" "${OUT_DIR}/nmf_k4_clinical_associations.csv"
echo "Wrote: nmf_k4_programs.csv, nmf_k4_program_labels.csv, nmf_k4_subtype_markers.csv, nmf_k4_subtype_pathways.csv, nmf_k4_clinical_associations.csv"

# Restore k=6 outputs from backup
echo "--- Restoring k=6 outputs ---"
for f in nmf_assignments.csv program_labels.csv subtype_markers.csv \
         subtype_pathways.csv subtype_clinical_associations.csv; do
  if [ -f "${OUT_DIR}/${f}.k6_backup" ]; then
    mv "${OUT_DIR}/${f}.k6_backup" "${OUT_DIR}/${f}"
  fi
done
echo "k=6 outputs restored."

echo ""
echo "=== Step 2 complete ==="

# ---- Step 3: k-means ARI sanity check ----
echo ""
echo "=== Step 3/3: Consensus k-means ARI sanity check ==="

export NMF_CACHE_PATH="RNA-seq/results/subtypes/nmf_results_cache_clean.rds"
Rscript RNA-seq/95b_nmf_k4_vs_kmeans_ari.R

echo ""
echo "=== Step 3 complete ==="

echo ""
echo "============================================================"
echo "  NMF k=4 co-primary pipeline COMPLETE"
echo "  End: $(date)"
echo ""
echo "  Output files:"
echo "    ${OUT_DIR}/nmf_k4_programs.csv"
echo "    ${OUT_DIR}/nmf_k4_seed_robustness.csv"
echo "    ${OUT_DIR}/nmf_vs_kmeans_ari.csv"
echo "    ${OUT_DIR}/nmf_k4_program_labels.csv"
echo "    ${OUT_DIR}/nmf_k4_subtype_markers.csv"
echo "    ${OUT_DIR}/nmf_k4_subtype_pathways.csv"
echo "    ${OUT_DIR}/nmf_k4_clinical_associations.csv"
echo "============================================================"
