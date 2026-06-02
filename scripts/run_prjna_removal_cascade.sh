#!/bin/bash
set -euo pipefail
# =====================================================================
# PRJNA512027 removal cascade — v2 (2026-05-15)
# ---------------------------------------------------------------------
# v1 failed at Step A with a sh-syntax error: --wrap runs through /bin/sh
# which can't eval micromamba's bash-only shell hook. v2 invokes
# `micromamba run -n rnaseq Rscript ...` directly via absolute path; no
# shell hook needed.
#
# Partition fan-out per user direction (cpu queue congested, io/bigmem/gpu
# idle):
#   * io      — light/short steps (A, B1, C4, D)
#   * cpu     — parallel mid-weight dream contrasts (B2-B6/B6b)
#   * bigmem  — atlas integrators (C1, C2, C3)
#   * gpu     — fallback only (cluster typically requires --gres=gpu)
#
# Dependency chain:
#   A (03_integrate)
#     ├─> B1 (05d_fibrosis)
#     ├─> B2 (13_nafl_vs_nash)
#     ├─> B3 (05e_mash_vs_masl_strict)
#     ├─> B4 (05f_mash_vs_healthy)
#     ├─> B5 (05g_masl_vs_healthy)
#     └─> B6 (26_sex_stratified) -> B6b (26b_sex_interaction)
#   {B1..B6b} -> C1 (27a_atlas) -> C2 (75_causal) -> C3 (217_strat_atlas)
#     -> C4 (48_cross_ancestry) -> D (run_pub_figures)
#
# Fallback policy:
#   * If QOSMaxMemoryPerUser blocks, resubmit step with --qos=interactive
#     (no memory cap, max 4 concurrent jobs).
#   * If io queue stalls, swap to --partition=cpu or --partition=gpu
#     --gres=gpu:1 --qos=interactive.
#   * If bigmem queue stalls for atlas, fall back --partition=cpu --mem=210G.
# =====================================================================

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
INT_DIR="${BASE}/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts"
RNA_DIR="${BASE}/RNA-seq"
FIG_RUN="${BASE}/scripts/figures/run_pub_figures.sh"
LOG_DIR="${BASE}/logs/prjna_removal_cascade"
mkdir -p "${LOG_DIR}"

MAMBA="/gpfs/commons/home/jameslee/.local/bin/micromamba"
EXPORTS="ALL,MASLD_PROJECT_ROOT=${BASE}"

# wrap-builder: produces a sh-safe --wrap payload
wrap_r() {
  local dir="$1" script="$2"
  printf 'set -eu; cd %q && %s run -n rnaseq Rscript %s\n' "$dir" "$MAMBA" "$script"
}

echo "=== PRJNA512027 removal cascade v2 — submission ==="
echo "Logs -> ${LOG_DIR}"
echo "Submitted at $(date)"
echo

# ---------------------------------------------------------------------
# Step A: 03_integrate_counts.R — io partition (idle)
# ---------------------------------------------------------------------
JOB_03=$(sbatch --parsable \
  --job-name=prjna_03_integrate \
  --partition=io --qos=nslab \
  --cpus-per-task=8 --mem=96G --time=48:00:00 \
  --export="${EXPORTS}" \
  --output="${LOG_DIR}/slurm-%j.out" --error="${LOG_DIR}/slurm-%j.out" \
  --wrap="$(wrap_r "${INT_DIR}" 03_integrate_counts.R)")
echo "Step A   [03_integrate_counts]         io      ${JOB_03}"

# ---------------------------------------------------------------------
# Step B1: 05d_fibrosis_vs_healthy_dream.R — io partition
# ---------------------------------------------------------------------
JOB_05D=$(sbatch --parsable \
  --dependency=afterok:${JOB_03} \
  --job-name=prjna_05d_fibrosis \
  --partition=io --qos=nslab \
  --cpus-per-task=8 --mem=64G --time=48:00:00 \
  --export="${EXPORTS}" \
  --output="${LOG_DIR}/slurm-%j.out" --error="${LOG_DIR}/slurm-%j.out" \
  --wrap="$(wrap_r "${INT_DIR}" 05d_fibrosis_vs_healthy_dream.R)")
echo "Step B1  [05d_fibrosis_vs_healthy]     io      ${JOB_05D}  (after ${JOB_03})"

# ---------------------------------------------------------------------
# Step B2: 13_nafl_vs_nash_de.R — cpu partition
# ---------------------------------------------------------------------
JOB_13=$(sbatch --parsable \
  --dependency=afterok:${JOB_03} \
  --job-name=prjna_13_nafl_vs_nash \
  --partition=cpu --qos=nslab \
  --cpus-per-task=8 --mem=96G --time=48:00:00 \
  --export="${EXPORTS}" \
  --output="${LOG_DIR}/slurm-%j.out" --error="${LOG_DIR}/slurm-%j.out" \
  --wrap="$(wrap_r "${INT_DIR}" 13_nafl_vs_nash_de.R)")
echo "Step B2  [13_nafl_vs_nash]             cpu     ${JOB_13}   (after ${JOB_03})"

# ---------------------------------------------------------------------
# Step B3: 05e_mash_vs_masl_dream_strict.R — cpu partition
# ---------------------------------------------------------------------
JOB_05E=$(sbatch --parsable \
  --dependency=afterok:${JOB_03} \
  --job-name=prjna_05e_mash_vs_masl \
  --partition=cpu --qos=nslab \
  --cpus-per-task=8 --mem=96G --time=48:00:00 \
  --export="${EXPORTS}" \
  --output="${LOG_DIR}/slurm-%j.out" --error="${LOG_DIR}/slurm-%j.out" \
  --wrap="$(wrap_r "${INT_DIR}" 05e_mash_vs_masl_dream_strict.R)")
echo "Step B3  [05e_mash_vs_masl_strict]     cpu     ${JOB_05E}  (after ${JOB_03})"

# ---------------------------------------------------------------------
# Step B4: 05f_mash_vs_healthy_dream.R — cpu partition
# ---------------------------------------------------------------------
JOB_05F=$(sbatch --parsable \
  --dependency=afterok:${JOB_03} \
  --job-name=prjna_05f_mash_vs_healthy \
  --partition=cpu --qos=nslab \
  --cpus-per-task=8 --mem=96G --time=48:00:00 \
  --export="${EXPORTS}" \
  --output="${LOG_DIR}/slurm-%j.out" --error="${LOG_DIR}/slurm-%j.out" \
  --wrap="$(wrap_r "${INT_DIR}" 05f_mash_vs_healthy_dream.R)")
echo "Step B4  [05f_mash_vs_healthy]         cpu     ${JOB_05F}  (after ${JOB_03})"

# ---------------------------------------------------------------------
# Step B5: 05g_masl_vs_healthy_dream.R — cpu partition
# ---------------------------------------------------------------------
JOB_05G=$(sbatch --parsable \
  --dependency=afterok:${JOB_03} \
  --job-name=prjna_05g_masl_vs_healthy \
  --partition=cpu --qos=nslab \
  --cpus-per-task=8 --mem=96G --time=48:00:00 \
  --export="${EXPORTS}" \
  --output="${LOG_DIR}/slurm-%j.out" --error="${LOG_DIR}/slurm-%j.out" \
  --wrap="$(wrap_r "${INT_DIR}" 05g_masl_vs_healthy_dream.R)")
echo "Step B5  [05g_masl_vs_healthy]         cpu     ${JOB_05G}  (after ${JOB_03})"

# ---------------------------------------------------------------------
# Step B6: 26_sex_stratified_analysis.R — cpu partition
# ---------------------------------------------------------------------
JOB_26=$(sbatch --parsable \
  --dependency=afterok:${JOB_03} \
  --job-name=prjna_26_sex_strat \
  --partition=cpu --qos=nslab \
  --cpus-per-task=8 --mem=96G --time=48:00:00 \
  --export="${EXPORTS}" \
  --output="${LOG_DIR}/slurm-%j.out" --error="${LOG_DIR}/slurm-%j.out" \
  --wrap="$(wrap_r "${INT_DIR}" 26_sex_stratified_analysis.R)")
echo "Step B6  [26_sex_stratified]           cpu     ${JOB_26}   (after ${JOB_03})"

# ---------------------------------------------------------------------
# Step B6b: 26b_sex_interaction_per_transition.R — cpu partition
# ---------------------------------------------------------------------
JOB_26B=$(sbatch --parsable \
  --dependency=afterok:${JOB_26} \
  --job-name=prjna_26b_sex_interaction \
  --partition=cpu --qos=nslab \
  --cpus-per-task=8 --mem=96G --time=48:00:00 \
  --export="${EXPORTS}" \
  --output="${LOG_DIR}/slurm-%j.out" --error="${LOG_DIR}/slurm-%j.out" \
  --wrap="$(wrap_r "${INT_DIR}" 26b_sex_interaction_per_transition.R)")
echo "Step B6b [26b_sex_interaction]         cpu     ${JOB_26B}  (after ${JOB_26})"

# ---------------------------------------------------------------------
# Step C1: 27a_assemble_evidence_atlas.R — bigmem
# ---------------------------------------------------------------------
JOB_27A=$(sbatch --parsable \
  --dependency=afterok:${JOB_05D}:${JOB_13}:${JOB_05E}:${JOB_05F}:${JOB_05G}:${JOB_26B} \
  --job-name=prjna_27a_atlas \
  --partition=bigmem --qos=nslab \
  --cpus-per-task=16 --mem=210G --time=72:00:00 \
  --export="${EXPORTS}" \
  --output="${LOG_DIR}/slurm-%j.out" --error="${LOG_DIR}/slurm-%j.out" \
  --wrap="$(wrap_r "${RNA_DIR}" 27a_assemble_evidence_atlas.R)")
echo "Step C1  [27a_assemble_atlas]          bigmem  ${JOB_27A}  (after B1-B6b)"

# ---------------------------------------------------------------------
# Step C2: 75_integrate_causal_overhaul.R — bigmem
# ---------------------------------------------------------------------
JOB_75=$(sbatch --parsable \
  --dependency=afterok:${JOB_27A} \
  --job-name=prjna_75_causal \
  --partition=bigmem --qos=nslab \
  --cpus-per-task=16 --mem=210G --time=72:00:00 \
  --export="${EXPORTS}" \
  --output="${LOG_DIR}/slurm-%j.out" --error="${LOG_DIR}/slurm-%j.out" \
  --wrap="$(wrap_r "${RNA_DIR}" 75_integrate_causal_overhaul.R)")
echo "Step C2  [75_integrate_causal]         bigmem  ${JOB_75}   (after ${JOB_27A})"

# ---------------------------------------------------------------------
# Step C3: 217_stratified_causal_atlas.R — bigmem
# ---------------------------------------------------------------------
JOB_217=$(sbatch --parsable \
  --dependency=afterok:${JOB_75} \
  --job-name=prjna_217_strat_atlas \
  --partition=bigmem --qos=nslab \
  --cpus-per-task=16 --mem=200G --time=72:00:00 \
  --export="${EXPORTS}" \
  --output="${LOG_DIR}/slurm-%j.out" --error="${LOG_DIR}/slurm-%j.out" \
  --wrap="$(wrap_r "${RNA_DIR}" 217_stratified_causal_atlas.R)")
echo "Step C3  [217_stratified_causal_atlas] bigmem  ${JOB_217}  (after ${JOB_75})"

# ---------------------------------------------------------------------
# Step C4: 48_cross_ancestry_replication.R — io partition
# ---------------------------------------------------------------------
JOB_48=$(sbatch --parsable \
  --dependency=afterok:${JOB_217} \
  --job-name=prjna_48_cross_anc \
  --partition=io --qos=nslab \
  --cpus-per-task=8 --mem=96G --time=48:00:00 \
  --export="${EXPORTS}" \
  --output="${LOG_DIR}/slurm-%j.out" --error="${LOG_DIR}/slurm-%j.out" \
  --wrap="$(wrap_r "${RNA_DIR}" 48_cross_ancestry_replication.R)")
echo "Step C4  [48_cross_ancestry]           io      ${JOB_48}   (after ${JOB_217})"

# ---------------------------------------------------------------------
# Step D: run_pub_figures.sh — io partition (overrides in-file headers)
# ---------------------------------------------------------------------
JOB_FIG=$(sbatch --parsable \
  --dependency=afterok:${JOB_48} \
  --job-name=prjna_pub_figures \
  --partition=io --qos=nslab \
  --cpus-per-task=4 --mem=64G --time=48:00:00 \
  --export="${EXPORTS}" \
  --output="${LOG_DIR}/slurm-%j.out" --error="${LOG_DIR}/slurm-%j.out" \
  "${FIG_RUN}")
echo "Step D   [run_pub_figures]             io      ${JOB_FIG}  (after ${JOB_48})"

echo
echo "=== All jobs submitted ==="
echo "Monitor with: squeue -u \$USER -o '%.10i %.30j %.8T %.10M %.6D %R'"
echo
echo "Job ID summary:"
echo "  A   03_integrate_counts          ${JOB_03}    io"
echo "  B1  05d_fibrosis_vs_healthy      ${JOB_05D}   io"
echo "  B2  13_nafl_vs_nash              ${JOB_13}    cpu"
echo "  B3  05e_mash_vs_masl_strict      ${JOB_05E}   cpu"
echo "  B4  05f_mash_vs_healthy          ${JOB_05F}   cpu"
echo "  B5  05g_masl_vs_healthy          ${JOB_05G}   cpu"
echo "  B6  26_sex_stratified            ${JOB_26}    cpu"
echo "  B6b 26b_sex_interaction          ${JOB_26B}   cpu"
echo "  C1  27a_assemble_atlas           ${JOB_27A}   bigmem"
echo "  C2  75_integrate_causal          ${JOB_75}    bigmem"
echo "  C3  217_stratified_causal_atlas  ${JOB_217}   bigmem"
echo "  C4  48_cross_ancestry            ${JOB_48}    io"
echo "  D   run_pub_figures              ${JOB_FIG}   io"
