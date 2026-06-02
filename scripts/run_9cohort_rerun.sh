#!/bin/bash
#SBATCH --job-name=liver_9cohort
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/9cohort_rerun_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/9cohort_rerun_%j.err
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=12:00:00

# 9-cohort re-run: Excludes PRJNA512027 (L0/S0 batch confounded with disease)
# and GSE167523 (no controls). Runs dream + all downstream scripts.
#
# Chain: 05 (dream) → 07 (consensus) → 08 (pathway) → 25 (deconv) → 26 (sex)
# Then: 27a+27b (atlas rebuild) submitted as dependent job by wrapper script.

set -euo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

export MASLD_PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "$MASLD_PROJECT_ROOT"
mkdir -p logs

INT_SCRIPTS="RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts"

echo "=========================================="
echo "  9-COHORT RE-RUN (excl PRJNA512027)"
echo "  Job: $SLURM_JOB_ID | CPUs: $SLURM_CPUS_PER_TASK"
echo "  Started: $(date)"
echo "=========================================="

# Script 05: Dream mega-analysis (9 cohorts)
echo ""
echo "=== 05: Dream Mega-Analysis (9 cohorts, 16 CPUs) ==="
Rscript "${INT_SCRIPTS}/05_dream_mega_analysis.R" 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 05_dream_mega_analysis.R"; exit 1; fi

# Script 07: Consensus DEGs
echo ""
echo "=== 07: Consensus DEGs ==="
Rscript "${INT_SCRIPTS}/07_consensus_degs.R" 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 07_consensus_degs.R"; exit 1; fi

# Script 08: Pathway analysis
echo ""
echo "=== 08: Pathway Analysis ==="
Rscript "${INT_SCRIPTS}/08_pathway_analysis.R" 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 08_pathway_analysis.R"; exit 1; fi

# Script 25: Deconvolution attribution (9 cohorts)
echo ""
echo "=== 25: Deconvolution Attribution ==="
Rscript "${INT_SCRIPTS}/25_deconv_attribution.R" 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 25_deconv_attribution.R"; exit 1; fi

# Script 26: Sex-stratified analysis (9 cohorts)
echo ""
echo "=== 26: Sex-Stratified Analysis ==="
Rscript "${INT_SCRIPTS}/26_sex_stratified_analysis.R" 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 26_sex_stratified_analysis.R"; exit 1; fi

echo ""
echo "=========================================="
echo "  Core pipeline complete: $(date)"
echo "=========================================="
