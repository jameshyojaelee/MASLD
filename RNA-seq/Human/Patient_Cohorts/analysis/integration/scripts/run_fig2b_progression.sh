#!/bin/bash
#SBATCH --job-name=fig2b_progression
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/fig2b_progression_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/fig2b_progression_%j.err
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=16
#SBATCH --mem=260G
#SBATCH --time=24:00:00

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPTS="${BASE}/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

export MASLD_PROJECT_ROOT="${BASE}"

echo "=== Fig 2b Progression Analyses: $(date) ==="
echo "SLURM_JOB_ID: ${SLURM_JOB_ID}"
echo "SLURM_CPUS_PER_TASK: ${SLURM_CPUS_PER_TASK}"
echo ""

# --- Analysis: Consecutive stage dream contrasts + landscape data ---
echo "=== 14d: Consecutive Stage Dream ==="
Rscript "${SCRIPTS}/14d_consecutive_stage_dream.R" 2>&1
if [ $? -ne 0 ]; then
  echo "FAILED: 14d_consecutive_stage_dream.R"
  exit 1
fi

echo ""

# --- Figure: Stage transitions (bar charts) ---
echo "=== Figure: Stage Transitions ==="
cd "${BASE}"
Rscript scripts/figures/figS_stage_transitions.R 2>&1
if [ $? -ne 0 ]; then
  echo "FAILED: figS_stage_transitions.R"
  exit 1
fi

echo ""

# --- Figure: Progression landscape (NAS × Fibrosis heatmap) ---
echo "=== Figure: Progression Landscape ==="
Rscript scripts/figures/figS_progression_landscape.R 2>&1
if [ $? -ne 0 ]; then
  echo "FAILED: figS_progression_landscape.R"
  exit 1
fi

echo ""
echo "=== All Fig 2b analyses completed: $(date) ==="
ls -lh "${BASE}/RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures/nas_consecutive_dream.csv" \
       "${BASE}/RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures/fibrosis_consecutive_dream.csv" \
       "${BASE}/RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures/progression_landscape_cells.csv" \
       2>/dev/null || true
ls -lh "${BASE}/figures/figS_stage_transitions.pdf" \
       "${BASE}/figures/figS_progression_landscape.pdf" \
       2>/dev/null || true
