#!/bin/bash
# run_nmf_k6_refit.sh — remediation R2: clean k=6 NMF RE-FIT (nrun=50) + relabel.
#
# Purpose: produce a fresh, post-C2-provenance k=6 NMF decomposition on the
# CURRENT expression matrix (merged_dge.rds, 2026-05-27) and a clean
# re-derivation of program labels / assignments / markers, so a single
# self-consistent k=6 labeling reconciles the 3-file discrepancy
# (program_labels.csv vs nmf_ksweep_program_atlas.csv vs manuscript).
#
# HARD GUARDRAIL: writes ONLY to NEW paths. The canonical cache + outputs
# (nmf_results_cache_clean.rds, program_labels.csv, nmf_assignments.csv,
# subtype_markers.csv) are NEVER touched.
#   - re-fit cache  : RNA-seq/results/subtypes/nmf_results_cache_clean_refit.rds
#   - shard         : RNA-seq/results/subtypes/nmf_shard_refit_k6.rds
#   - relabel outputs: RNA-seq/results/subtypes_refit/{program_labels,nmf_assignments,subtype_markers,...}.csv
#
# Matrix-build recipe MIRRORS the canonical clean cache (Script 95):
#   NMF_STRIP_HLA=0, NMF_PROTEIN_CODING_ONLY=0, N_TOP_GENES=5000, top-IQR.
# k=6 fit MIRRORS canonical: nrun=50, brunet, seed=random, set.seed(42+6).
#
# Chain (sbatch --dependency=afterok): 95 matrix -> 94a k=6 -> 94b merge -> 44 relabel.
# Job-name: NMF (single word; user tracks via squeue --name=NMF).
# Partition: io + --qos=interactive (cpu is full + nslab QOS mem-capped 2026-06-10).
#
# Submit: bash RNA-seq/run_nmf_k6_refit.sh
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

REFIT_CACHE="RNA-seq/results/subtypes/nmf_results_cache_clean_refit.rds"
SHARD_TAG="refit"
REFIT_OUT="RNA-seq/results/subtypes_refit"
LOG_DIR="RNA-seq/logs"
TMP_DIR="RNA-seq/logs/nmf_k6_refit_scripts"
mkdir -p "$LOG_DIR" "$TMP_DIR" "$REFIT_OUT"

INIT='eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)" && micromamba activate rnaseq && cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design'

# ── Step 1: Build clean matrix into the NEW refit cache (Script 95) ──────────
# Script 95 sees no cache at REFIT_CACHE -> rebuilds mat_nn/mat from
# merged_dge.rds with the canonical strip, writes cache shell (+ inline k=3,4),
# then exits cleanly (k=5..8 rubric skipped — we only need the matrix here).
cat > "$TMP_DIR/step1_matrix.sh" << EOF
#!/bin/bash
#SBATCH --partition=cpu --cpus-per-task=8 --mem=80G --time=48:00:00
#SBATCH --job-name=NMF
#SBATCH --output=$LOG_DIR/nmf_k6refit_matrix_%j.out
#SBATCH --error=$LOG_DIR/nmf_k6refit_matrix_%j.err
set -e
$INIT
export NMF_CACHE_PATH="$REFIT_CACHE"
export NMF_SHARD_TAG="$SHARD_TAG"
export NMF_STRIP_HLA=0
export NMF_PROTEIN_CODING_ONLY=0
export NMF_RUNS=50
Rscript RNA-seq/95_nmf_clean_ksweep.R
EOF

echo "Submitting Step 1: matrix build (io/interactive) ..."
JOB1=$(sbatch --parsable "$TMP_DIR/step1_matrix.sh")
echo "  matrix job: $JOB1"

# ── Step 2: Fit k=6 (nrun=50) as a refit shard (Script 94a) ──────────────────
cat > "$TMP_DIR/step2_k6.sh" << EOF
#!/bin/bash
#SBATCH --partition=cpu --cpus-per-task=8 --mem=64G --time=48:00:00
#SBATCH --job-name=NMF
#SBATCH --output=$LOG_DIR/nmf_k6refit_k6_%j.out
#SBATCH --error=$LOG_DIR/nmf_k6refit_k6_%j.err
set -e
$INIT
export NMF_K=6
export NMF_RUNS=50
export NMF_CACHE_PATH="$REFIT_CACHE"
export NMF_SHARD_TAG="$SHARD_TAG"
Rscript RNA-seq/94a_nmf_single_k.R
EOF

echo "Submitting Step 2: k=6 fit (nrun=50) ..."
JOB2=$(sbatch --parsable --dependency=afterok:"$JOB1" "$TMP_DIR/step2_k6.sh")
echo "  k=6 fit job: $JOB2"

# ── Step 3: Merge k=6 shard into the refit cache (Script 94b) ────────────────
cat > "$TMP_DIR/step3_merge.sh" << EOF
#!/bin/bash
#SBATCH --partition=cpu --cpus-per-task=2 --mem=32G --time=48:00:00
#SBATCH --job-name=NMF
#SBATCH --output=$LOG_DIR/nmf_k6refit_merge_%j.out
#SBATCH --error=$LOG_DIR/nmf_k6refit_merge_%j.err
set -e
$INIT
export NMF_CACHE_PATH="$REFIT_CACHE"
export NMF_SHARD_TAG="$SHARD_TAG"
Rscript RNA-seq/94b_merge_shards.R
EOF

echo "Submitting Step 3: merge ..."
JOB3=$(sbatch --parsable --dependency=afterok:"$JOB2" "$TMP_DIR/step3_merge.sh")
echo "  merge job: $JOB3"

# ── Step 4: Relabel k=6 (Script 44) -> subtypes_refit/ (canonical untouched) ──
cat > "$TMP_DIR/step4_relabel.sh" << EOF
#!/bin/bash
#SBATCH --partition=cpu --cpus-per-task=8 --mem=80G --time=48:00:00
#SBATCH --job-name=NMF
#SBATCH --output=$LOG_DIR/nmf_k6refit_relabel_%j.out
#SBATCH --error=$LOG_DIR/nmf_k6refit_relabel_%j.err
set -e
$INIT
export NMF_CHOSEN_K=6
export NMF_CACHE_PATH="$REFIT_CACHE"
export SUBTYPE_OUT_DIR="$REFIT_OUT"
Rscript RNA-seq/44_molecular_subtyping.R
EOF

echo "Submitting Step 4: relabel ..."
JOB4=$(sbatch --parsable --dependency=afterok:"$JOB3" "$TMP_DIR/step4_relabel.sh")
echo "  relabel job: $JOB4"

echo ""
echo "Chain: matrix($JOB1) -> k6_fit($JOB2) -> merge($JOB3) -> relabel($JOB4)"
echo "Refit cache : $REFIT_CACHE"
echo "Refit outputs: $REFIT_OUT/"
echo "Monitor: squeue -u \$USER --name=NMF"
