#!/bin/bash
#SBATCH --job-name=fig2_age_sex
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/fig2_age_sex_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/fig2_age_sex_%j.err
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=210G
#SBATCH --time=12:00:00

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPTS="${BASE}/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts"
LOG_DIR="${BASE}/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs"

mkdir -p "${LOG_DIR}"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

export MASLD_PROJECT_ROOT="${BASE}"

echo "=== Fig 2 Age+Sex Analyses: $(date) ==="
echo "SLURM_JOB_ID: ${SLURM_JOB_ID}"
echo "SLURM_CPUS_PER_TASK: ${SLURM_CPUS_PER_TASK}"
echo ""

# --- Step 1: Age + sex stratified dream (bigmem, ~4h) ---
echo "=== 14e: Age + Sex Stratified Dream ==="
Rscript "${SCRIPTS}/14e_age_sex_stratified_dream.R" 2>&1
if [ $? -ne 0 ]; then
  echo "FAILED: 14e_age_sex_stratified_dream.R"
  exit 1
fi
echo ""

# --- Step 2: NAS score GSEA (lightweight, ~30min) ---
echo "=== 14f: NAS Score GSEA ==="
Rscript "${SCRIPTS}/14f_nas_stage_gsea.R" 2>&1
if [ $? -ne 0 ]; then
  echo "FAILED: 14f_nas_stage_gsea.R"
  exit 1
fi
echo ""

# --- Step 3: Generate Figure 2 ---
echo "=== Figure 2: Disease Progression ==="
cd "${BASE}"
Rscript scripts/figures/fig2_disease_progression.R 2>&1
if [ $? -ne 0 ]; then
  echo "FAILED: fig2_disease_progression.R"
  exit 1
fi
echo ""

# --- Step 4: Generate Supplementary Figure 8 ---
echo "=== Supplementary Figure 8 ==="
Rscript scripts/figures/figS8_disease_progression_supp.R 2>&1
if [ $? -ne 0 ]; then
  echo "FAILED: figS8_disease_progression_supp.R"
  exit 1
fi
echo ""

echo "=== All Fig 2 analyses completed: $(date) ==="
echo ""
echo "Output files:"
ls -lh "${BASE}/RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures/per_sample_dysregulation.csv" \
       "${BASE}/RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures/age_stratified_dream.csv" \
       "${BASE}/RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures/nas_score_gsea.csv" \
       "${BASE}/figures/fig2_disease_progression.pdf" \
       "${BASE}/figures/figS8_disease_progression_supp.pdf" \
       2>/dev/null || true
