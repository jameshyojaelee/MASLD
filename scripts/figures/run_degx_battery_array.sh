#!/bin/bash
#SBATCH --job-name=Rexact
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/rexact_shard_%A_%a.log
#SBATCH --array=0-18
# Parallel fresh re-run of the full degx R-exact battery: one method per array task,
# each into its own shard dir, so the whole battery finishes in ~slowest-method time
# instead of ~3h serial. Gather + compare-vs-stored-degx happens after.

set -euo pipefail

PROJECT_ROOT="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
DEGX_HOME=/gpfs/commons/home/jameslee/degx

# 19 disease_vs_control method variants (exactly those in the stored degx manifest)
METHODS=(deseq2_wald deseq2_lrt deseq2_apeglm deseq2_ashr \
         edger_qlf edger_qlf_robust edger_exact \
         limma_voom limma_trend limma_voom_qw dream \
         metafor_deseq2_re metafor_deseq2_fe metafor_deseq2_hk \
         metafor_voom_re metafor_voom_fe metafor_voom_hk \
         combatseq_deseq2 sva_limma)

M="${METHODS[$SLURM_ARRAY_TASK_ID]}"
# Shards go to a work dir so deg_tables/ only ever holds the final consistent
# deg_<method>_degx.csv files (compare_full_battery_vs_degx.R gathers from here).
OUT="$PROJECT_ROOT/figures/supplementary/figS_methods_validation/disease_signature_sweep/work/degx_rerun/shards/$M"
mkdir -p "$PROJECT_ROOT/logs" "$OUT"
cd "$DEGX_HOME"   # engine sources R/corrections.R relative to getwd()

echo "[$(date)] shard ${SLURM_ARRAY_TASK_ID} method=$M on $(hostname)"
micromamba run -n rnaseq Rscript R/run_methods_real.R \
  --contrast disease_vs_control \
  --data_dir "$DEGX_HOME/data/masld" \
  --out_dir  "$OUT" \
  --methods  "$M"
echo "[$(date)] done $M -> $OUT"
