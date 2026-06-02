#!/bin/bash
# run_nmf_v2_pipeline.sh
# Full NMF v2 rerun on a fresh cache with the pseudogene strip applied.
# v1 cache and outputs left intact as backup.
# Partitions: io (matrix/merge/commit), gpu (compute shards — cpu is busy).
#
# Submit: bash RNA-seq/run_nmf_v2_pipeline.sh
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

V2_CACHE="RNA-seq/results/subtypes/nmf_results_cache_clean_v2.rds"
SHARD_TAG="clean_v2"
LOG_DIR="RNA-seq/logs"
TMP_DIR="RNA-seq/logs/nmf_v2_scripts"
mkdir -p "$LOG_DIR" "$TMP_DIR"

INIT='eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)" && micromamba activate rnaseq && cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design'

# ── Step 1: Build matrix + k=3,4 (io — cpu is busy) ──────────────────────
cat > "$TMP_DIR/step1_matrix.sh" << EOF
#!/bin/bash
#SBATCH --partition=io --qos=interactive --cpus-per-task=8 --mem=80G --time=72:00:00
#SBATCH --job-name=nmf_v2_matrix
#SBATCH --output=$LOG_DIR/nmf_v2_matrix_%j.out
#SBATCH --error=$LOG_DIR/nmf_v2_matrix_%j.err
set -e
$INIT
export NMF_CACHE_PATH="$V2_CACHE"
export NMF_SHARD_TAG="$SHARD_TAG"
export NMF_STRIP_HLA=0
export NMF_PROTEIN_CODING_ONLY=0
Rscript RNA-seq/95_nmf_clean_ksweep.R
EOF

echo "Submitting Step 1: matrix build ..."
JOB1=$(sbatch --parsable "$TMP_DIR/step1_matrix.sh")
echo "  job: $JOB1"

# ── Step 2: Shards k=5..8 in parallel (gpu — idle CPUs available) ─────────
echo "Submitting Step 2: shards k=5..8 (gpu) ..."
SHARD_JOBS=()
for K in 5 6 7 8; do
  cat > "$TMP_DIR/step2_shard_k${K}.sh" << EOF
#!/bin/bash
#SBATCH --partition=gpu --cpus-per-task=8 --mem=64G --time=48:00:00
#SBATCH --job-name=nmf_v2_k${K}
#SBATCH --output=$LOG_DIR/nmf_v2_k${K}_%j.out
#SBATCH --error=$LOG_DIR/nmf_v2_k${K}_%j.err
set -e
$INIT
export NMF_K=$K
export NMF_CACHE_PATH="$V2_CACHE"
export NMF_SHARD_TAG="$SHARD_TAG"
Rscript RNA-seq/94a_nmf_single_k.R
EOF
  JID=$(sbatch --parsable --dependency=afterok:"$JOB1" "$TMP_DIR/step2_shard_k${K}.sh")
  SHARD_JOBS+=("$JID")
  echo "  k=$K job: $JID"
done

SHARD_DEP="afterok:$(IFS=:; echo "${SHARD_JOBS[*]}")"

# ── Step 3: Merge shards (io — lightweight) ───────────────────────────────
cat > "$TMP_DIR/step3_merge.sh" << EOF
#!/bin/bash
#SBATCH --partition=io --qos=interactive --cpus-per-task=2 --mem=32G --time=4:00:00
#SBATCH --job-name=nmf_v2_merge
#SBATCH --output=$LOG_DIR/nmf_v2_merge_%j.out
#SBATCH --error=$LOG_DIR/nmf_v2_merge_%j.err
set -e
$INIT
export NMF_CACHE_PATH="$V2_CACHE"
export NMF_SHARD_TAG="$SHARD_TAG"
Rscript RNA-seq/94b_merge_shards.R
EOF

echo "Submitting Step 3: merge ..."
JOB3=$(sbatch --parsable --dependency="$SHARD_DEP" "$TMP_DIR/step3_merge.sh")
echo "  job: $JOB3"

# ── Step 4: Rubric + relabel + downstream + figures (io) ──────────────────
cat > "$TMP_DIR/step4_commit.sh" << EOF
#!/bin/bash
#SBATCH --partition=io --qos=interactive --cpus-per-task=8 --mem=80G --time=72:00:00
#SBATCH --job-name=nmf_v2_commit
#SBATCH --output=$LOG_DIR/nmf_v2_commit_%j.out
#SBATCH --error=$LOG_DIR/nmf_v2_commit_%j.err
set -e
$INIT
export NMF_CACHE_PATH="$V2_CACHE"
export NMF_SHARD_TAG="$SHARD_TAG"
export NMF_STRIP_HLA=0
export NMF_PROTEIN_CODING_ONLY=0

echo "=== 4a: rubric ==="
Rscript RNA-seq/95_nmf_clean_ksweep.R

echo "=== 4b: relabel ==="
Rscript RNA-seq/95c_relabel_programs.R

echo "=== 4c: assignments + atlas ==="
Rscript RNA-seq/44_molecular_subtyping.R
Rscript RNA-seq/208_subtype_coloc.R
Rscript RNA-seq/217_stratified_causal_atlas.R

echo "=== 4d: figures ==="
Rscript scripts/figures/fig2_panel_nmf_programs.R
Rscript scripts/figures/fig2_panel_nmf_programs_lines.R
Rscript scripts/figures/figS08_nmf_subtyping.R

echo "=== NMF v2 complete ==="
EOF

echo "Submitting Step 4: commit ..."
JOB4=$(sbatch --parsable --dependency=afterok:"$JOB3" "$TMP_DIR/step4_commit.sh")
echo "  job: $JOB4"

echo ""
echo "Chain: matrix($JOB1) → shards(${SHARD_JOBS[*]}) → merge($JOB3) → commit($JOB4)"
echo "Monitor: squeue -u \$USER --name=nmf_v2_matrix,nmf_v2_k5,nmf_v2_k6,nmf_v2_k7,nmf_v2_k8,nmf_v2_merge,nmf_v2_commit"
