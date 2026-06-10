#!/usr/bin/env bash
# ============================================================================
# WS3_01_refit_lmm_augmented.sh
#
# Re-fit the stage-LR mixed model (Script 346) on the F_stage_augmented axis
# (fine axis for WS3 stage-DEG routing). The canonical gated fstage LMM
# (stage_lr_lmm_fstage.tsv) was fit on F_stage_inferred because the
# BOOTSTRAP_FALLBACK_REQUIRED.flag is present. WS3 explicitly wants the
# augmented fine axis, so we force ALLOW_F_STAGE_AUGMENTED_DESPITE_FLAG=TRUE.
#
# NON-DESTRUCTIVE: backs up the canonical fstage LMM before the run, then
# moves the augmented output to stage_lr_lmm_fstage_augmented.tsv and restores
# the canonical inferred-axis file. The coarse output is identical re-computed
# (Healthy-reference contrasts) so we leave it. CPU only, reuses parquets.
#
# Usage:
#   env -u SLURM_JOB_ID sbatch --partition=cpu --qos=interactive --mem=48G \
#       --time=48:00:00 --job-name=liana \
#       Analysis/SingleCell/scripts/WS3_01_refit_lmm_augmented.sh
# ============================================================================
#SBATCH --output=Analysis/SingleCell/scripts/logs/WS3_01_refit_lmm_augmented_%j.out
#SBATCH --error=Analysis/SingleCell/scripts/logs/WS3_01_refit_lmm_augmented_%j.err
#SBATCH --cpus-per-task=16

set -eo pipefail

BASE="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
cd "$BASE"
mkdir -p Analysis/SingleCell/scripts/logs

D="$BASE/Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory"
CANON="$D/stage_lr_lmm_fstage.tsv"
CANON_BAK="$D/stage_lr_lmm_fstage.tsv.canon_inferred_bak"
AUG_OUT="$D/stage_lr_lmm_fstage_augmented.tsv"
COARSE_BAK="$D/stage_lr_lmm_coarse.tsv.bak_ws3"

echo "[WS3] activating rnaseq env"
export PATH="/gpfs/commons/home/jameslee/.local/bin:$PATH"
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u   # enable unbound-var check only after env activate.d scripts have run

# Preserve canonical gated outputs so the override re-run does not clobber them.
if [[ -f "$CANON" && ! -f "$CANON_BAK" ]]; then
  cp -p "$CANON" "$CANON_BAK"
  echo "[WS3] backed up canonical inferred-axis fstage LMM -> $CANON_BAK"
fi
if [[ -f "$D/stage_lr_lmm_coarse.tsv" && ! -f "$COARSE_BAK" ]]; then
  cp -p "$D/stage_lr_lmm_coarse.tsv" "$COARSE_BAK"
fi

echo "[WS3] re-running Script 346 with ALLOW_F_STAGE_AUGMENTED_DESPITE_FLAG=TRUE"
export ALLOW_F_STAGE_AUGMENTED_DESPITE_FLAG=TRUE
export SLURM_CPUS_PER_TASK="${SLURM_CPUS_PER_TASK:-16}"

Rscript "$BASE/Analysis/SingleCell/scripts/346_stage_ccc_mixedmodel.R"

# The override run wrote the augmented axis to stage_lr_lmm_fstage.tsv.
# Move it to the augmented name, then restore the canonical inferred file.
if [[ -f "$CANON" ]]; then
  mv "$CANON" "$AUG_OUT"
  echo "[WS3] augmented-axis fstage LMM -> $AUG_OUT"
fi
if [[ -f "$CANON_BAK" ]]; then
  cp -p "$CANON_BAK" "$CANON"
  echo "[WS3] restored canonical inferred-axis fstage LMM -> $CANON"
fi
# Restore canonical coarse (the re-run reproduces it identically, but be safe).
if [[ -f "$COARSE_BAK" ]]; then
  cp -p "$COARSE_BAK" "$D/stage_lr_lmm_coarse.tsv"
fi

echo "[WS3] augmented LMM re-fit complete."
echo "[WS3] augmented fstage rows: $(wc -l < "$AUG_OUT")"
