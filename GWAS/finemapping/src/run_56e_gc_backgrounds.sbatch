#!/bin/bash
#SBATCH --job-name=mb_gc_bg
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --array=0-2
#SBATCH --output=../logs/56e_bg_%A_%a.out
#SBATCH --error=../logs/56e_bg_%A_%a.err

# A4: GC-aware nucleotide backgrounds for motifbreakR
# Array index 0,1,2 -> uniform / genome / peak.
# Task 0 computes the BG files (under a small lock so tasks 1/2 don't race).
# Then each task runs Script 56_v2 with MB_BACKGROUND=<bg>.
#
# Submit:
#   sbatch run_56e_gc_backgrounds.sbatch
# Then submit the comparison as a dependent job (see footer).

set -eo pipefail

BG_MODES=("uniform" "genome" "peak")
export MB_BACKGROUND="${BG_MODES[$SLURM_ARRAY_TASK_ID]}"

echo "============================================================"
echo "Array task $SLURM_ARRAY_TASK_ID  | MB_BACKGROUND=$MB_BACKGROUND"
echo "Host: $(hostname)  | Date: $(date)"
echo "============================================================"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate motifbreakr
module load meme/5.5.9 || echo "meme module load skipped (FIMO disabled in v2)"

SRC_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/src"
RES_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/gwas_atac"
LOG_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs"
cd "$SRC_DIR"
mkdir -p "$LOG_DIR" "$RES_DIR"

# Background files needed by genome/peak tasks. Task 0 (uniform) computes them
# first and writes a flag file; tasks 1 and 2 wait for the flag.
FLAG="$RES_DIR/.bg_ready"
GENOME_BG="$RES_DIR/motifbreakr_bg_genome.rds"
PEAK_BG="$RES_DIR/motifbreakr_bg_peaks.rds"

if [[ "$SLURM_ARRAY_TASK_ID" == "0" ]]; then
  if [[ -f "$GENOME_BG" && -f "$PEAK_BG" ]]; then
    echo "[task 0] background files already exist — skipping 56e_gc_background.R"
  else
    echo "[task 0] computing nucleotide backgrounds (56e_gc_background.R) ..."
    Rscript 56e_gc_background.R
  fi
  touch "$FLAG"
else
  echo "[task $SLURM_ARRAY_TASK_ID] waiting for $FLAG ..."
  for i in $(seq 1 600); do
    if [[ -f "$FLAG" && -f "$GENOME_BG" && -f "$PEAK_BG" ]]; then
      echo "[task $SLURM_ARRAY_TASK_ID] backgrounds ready after ${i}0s"
      break
    fi
    sleep 10
  done
  if [[ ! -f "$GENOME_BG" || ! -f "$PEAK_BG" ]]; then
    echo "ERROR: backgrounds not produced within 100 min."
    exit 1
  fi
fi

echo "Running 56_motif_disruption_v2.R with bg=$MB_BACKGROUND ..."
Rscript 56_motif_disruption_v2.R

echo "Task $SLURM_ARRAY_TASK_ID done at $(date)"

# After array completes, submit the comparison job:
#   sbatch --dependency=afterok:$ARRAY_JOB_ID --partition=cpu --mem=32G \
#          --time=12:00:00 --job-name=mb_bg_cmp \
#          --output=../logs/56e_cmp_%j.out --error=../logs/56e_cmp_%j.err \
#          --wrap="source ~/.bashrc && micromamba activate rnaseq && \
#                  cd $SRC_DIR && Rscript 56e_compare_backgrounds.R"
