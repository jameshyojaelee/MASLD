#!/bin/bash
#SBATCH --job-name=fig2_disease_prog
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/fig2_analyses_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/fig2_analyses_%j.err
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

echo "=== Fig 2 Analyses: $(date) ==="
echo "SLURM_JOB_ID: ${SLURM_JOB_ID}"
echo "SLURM_CPUS_PER_TASK: ${SLURM_CPUS_PER_TASK}"
echo ""

# --- Analysis 1+2: Per-NAS-score + Per-fibrosis-stage dream ---
echo "=== 14b: NAS + Fibrosis Stage Dream ==="
Rscript "${SCRIPTS}/14b_nas_fibrosis_stage_dream.R" 2>&1
if [ $? -ne 0 ]; then
  echo "FAILED: 14b_nas_fibrosis_stage_dream.R"
  exit 1
fi

echo ""

# --- Analysis 3: Per-stage GSEA ---
echo "=== 14c: Per-Stage GSEA ==="
Rscript "${SCRIPTS}/14c_fibrosis_stage_gsea.R" 2>&1
if [ $? -ne 0 ]; then
  echo "FAILED: 14c_fibrosis_stage_gsea.R"
  exit 1
fi

echo ""

# --- Generate Figure 2 ---
echo "=== Figure 2: Disease Progression ==="
cd "${BASE}"
Rscript scripts/figures/fig2_disease_progression.R 2>&1
if [ $? -ne 0 ]; then
  echo "FAILED: fig2_disease_progression.R"
  exit 1
fi

echo ""
echo "=== All Fig 2 analyses completed: $(date) ==="
ls -lh "${BASE}/RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures/nas_score_dream.csv" \
       "${BASE}/RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures/fibrosis_stage_dream.csv" \
       "${BASE}/RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures/fibrosis_stage_gsea.csv" \
       2>/dev/null || true
ls -lh "${BASE}/figures/fig2_disease_progression.pdf" 2>/dev/null || true
