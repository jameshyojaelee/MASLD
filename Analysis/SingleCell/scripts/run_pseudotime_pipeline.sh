#!/bin/bash
set -euo pipefail

# =============================================================================
# Pseudotime Pipeline Orchestrator (Scripts 300-307)
#
# DAG:
#   P1 (300_preprocessing.py)
#     -> P2 (301_palantir_dpt.py)        \
#     -> P3 (302_cytotrace2.py)           -> P7 (306_consensus.py) -> P8 (307_figures.R)
#     -> P4 (303_paga_cross_celltype.py)  |
#     -> P5 (304_monocle3_trajectory.R)   -> P6 (305_bulk_sc.py) -/
#
# P2, P3, P4, P5 run in parallel after P1.
# P6 depends on P2 AND P5.
# P7 depends on P2 AND P3 AND P5.
# P8 depends on P7.
# =============================================================================

export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
BASE="$MASLD_PROJECT_ROOT"
SC_DIR="$BASE/Analysis/SingleCell"
SCRIPT_DIR="$SC_DIR/scripts"
LOGDIR="$SCRIPT_DIR/logs/pseudotime"

mkdir -p "$LOGDIR"

echo "=== Pseudotime Pipeline ==="
echo "Start: $(date)"
echo "Logs:  $LOGDIR"
echo ""

# --- Shared environment preambles ---
GPU_PREAMBLE='
set -euo pipefail
eval "$(micromamba shell hook --shell=bash)"
micromamba activate rapids_singlecell
export LD_LIBRARY_PATH=$(python3 -c "import nvidia.cusparselt.lib; print(nvidia.cusparselt.lib.__path__[0])"):${LD_LIBRARY_PATH:-}
export MASLD_PROJECT_ROOT='"$BASE"'
cd '"$BASE"'
'

CPU_PY_PREAMBLE='
set -euo pipefail
eval "$(micromamba shell hook --shell=bash)"
micromamba activate rapids_singlecell
export MASLD_PROJECT_ROOT='"$BASE"'
cd '"$BASE"'
'

R_PREAMBLE='
set -euo pipefail
eval "$(micromamba shell hook --shell=bash)"
micromamba activate rnaseq
export MASLD_PROJECT_ROOT='"$BASE"'
cd '"$BASE"'
'

# =============================================================================
# P1: Preprocessing (bigmem — 1.23M cell atlas needs >128GB)
# =============================================================================
JOB1=$(sbatch --parsable \
  --job-name=pt_300 \
  --partition=bigmem \
  --cpus-per-task=16 \
  --mem=256G \
  --time=48:00:00 \
  --output="$LOGDIR/pt_300_%j.log" \
  --wrap="bash -c '${CPU_PY_PREAMBLE}
python3 ${SCRIPT_DIR}/300_pseudotime_preprocessing.py
'")
echo "P1 [300_preprocessing]:       SLURM $JOB1"

# =============================================================================
# P2: Palantir DPT (GPU) — depends on P1
# NOTE: palantir + statsmodels must be pre-installed in rapids_singlecell env
# =============================================================================
JOB2=$(sbatch --parsable \
  --dependency=afterok:$JOB1 \
  --job-name=pt_301 \
  --partition=gpu \
  --gres=gpu:1 \
  --cpus-per-task=16 \
  --mem=128G \
  --time=48:00:00 \
  --output="$LOGDIR/pt_301_%j.log" \
  --wrap="bash -c '${GPU_PREAMBLE}
python3 ${SCRIPT_DIR}/301_palantir_dpt.py
'")
echo "P2 [301_palantir_dpt]:        SLURM $JOB2 (after P1: $JOB1)"

# =============================================================================
# P3: CytoTRACE 2 (GPU) — depends on P1 only
# =============================================================================
JOB3=$(sbatch --parsable \
  --dependency=afterok:$JOB1 \
  --job-name=pt_302 \
  --partition=gpu \
  --gres=gpu:1 \
  --cpus-per-task=16 \
  --mem=128G \
  --time=48:00:00 \
  --output="$LOGDIR/pt_302_%j.log" \
  --wrap="bash -c '${GPU_PREAMBLE}
python3 ${SCRIPT_DIR}/302_cytotrace2.py
'")
echo "P3 [302_cytotrace2]:          SLURM $JOB3 (after P1: $JOB1)"

# =============================================================================
# P4: PAGA cross-celltype (bigmem — loads full atlas)
# =============================================================================
JOB4=$(sbatch --parsable \
  --dependency=afterok:$JOB1 \
  --job-name=pt_303 \
  --partition=bigmem \
  --cpus-per-task=16 \
  --mem=256G \
  --time=48:00:00 \
  --output="$LOGDIR/pt_303_%j.log" \
  --wrap="bash -c '${CPU_PY_PREAMBLE}
python3 ${SCRIPT_DIR}/303_paga_cross_celltype.py
'")
echo "P4 [303_paga_cross_celltype]: SLURM $JOB4 (after P1: $JOB1)"

# =============================================================================
# P5: Monocle3 trajectory (R/CPU) — depends on P1 only
# =============================================================================
JOB5=$(sbatch --parsable \
  --dependency=afterok:$JOB1 \
  --job-name=pt_304 \
  --partition=cpu \
  --cpus-per-task=16 \
  --mem=128G \
  --time=48:00:00 \
  --output="$LOGDIR/pt_304_%j.log" \
  --wrap="bash -c '${R_PREAMBLE}
Rscript ${SCRIPT_DIR}/304_monocle3_trajectory.R
'")
echo "P5 [304_monocle3_trajectory]: SLURM $JOB5 (after P1: $JOB1)"
echo ""

# =============================================================================
# P6: Bulk-SC integration (GPU) — depends on P2 AND P5
# =============================================================================
JOB6=$(sbatch --parsable \
  --dependency=afterok:$JOB2:$JOB5 \
  --job-name=pt_305 \
  --partition=gpu \
  --gres=gpu:1 \
  --cpus-per-task=16 \
  --mem=128G \
  --time=48:00:00 \
  --output="$LOGDIR/pt_305_%j.log" \
  --wrap="bash -c '${GPU_PREAMBLE}
python3 ${SCRIPT_DIR}/305_bulk_sc_integration.py
'")
echo "P6 [305_bulk_sc]:             SLURM $JOB6 (after P2: $JOB2, P5: $JOB5)"

# =============================================================================
# P7: Consensus (GPU) — depends on P2 AND P3 AND P5
# =============================================================================
JOB7=$(sbatch --parsable \
  --dependency=afterok:$JOB2:$JOB3:$JOB5 \
  --job-name=pt_306 \
  --partition=gpu \
  --gres=gpu:1 \
  --cpus-per-task=16 \
  --mem=128G \
  --time=48:00:00 \
  --output="$LOGDIR/pt_306_%j.log" \
  --wrap="bash -c '${GPU_PREAMBLE}
python3 ${SCRIPT_DIR}/306_consensus_pseudotime.py
'")
echo "P7 [306_consensus]:           SLURM $JOB7 (after P2: $JOB2, P3: $JOB3, P5: $JOB5)"

# =============================================================================
# P8: Figures (R/CPU) — depends on P7
# =============================================================================
JOB8=$(sbatch --parsable \
  --dependency=afterok:$JOB7 \
  --job-name=pt_307 \
  --partition=cpu \
  --cpus-per-task=16 \
  --mem=128G \
  --time=48:00:00 \
  --output="$LOGDIR/pt_307_%j.log" \
  --wrap="bash -c '${R_PREAMBLE}
Rscript ${SCRIPT_DIR}/307_pseudotime_figures.R
'")
echo "P8 [307_figures]:             SLURM $JOB8 (after P7: $JOB7)"

# =============================================================================
# Summary
# =============================================================================
echo ""
echo "=== Pipeline submitted ==="
echo "P1  [300_preprocessing]:       $JOB1"
echo "P2  [301_palantir_dpt]:        $JOB2"
echo "P3  [302_cytotrace2]:          $JOB3"
echo "P4  [303_paga_cross_celltype]: $JOB4"
echo "P5  [304_monocle3_trajectory]: $JOB5"
echo "P6  [305_bulk_sc]:             $JOB6"
echo "P7  [306_consensus]:           $JOB7"
echo "P8  [307_figures]:             $JOB8"
echo ""
echo "Monitor: squeue -u $USER -j $JOB1,$JOB2,$JOB3,$JOB4,$JOB5,$JOB6,$JOB7,$JOB8"
echo "Logs:    $LOGDIR"
