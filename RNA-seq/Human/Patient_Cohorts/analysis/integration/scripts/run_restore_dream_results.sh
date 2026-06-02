#!/bin/bash
#SBATCH --job-name=liver_restore_dream05
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=32
#SBATCH --mem=300G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/restore_dream05_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/restore_dream05_%j.err

# PURPOSE: Restore dream_results.csv after dream_variants "base" run
# (which used OLD broken merged_dge.rds) overwrote the correct 18,774 DEG result.
# merged_dge.rds is correct (written 06:12 by Script 03 fix).
# consensus_degs.csv is correct (written 07:21 based on 18,774 DEGs).
# This job re-runs Script 05 on the correct data to restore dream_results.csv,
# then confirms consensus_degs.csv is still correct.

set -euo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

SCRIPTS_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts"
cd "$SCRIPTS_DIR"

echo "=== Restore dream_results.csv ==="
echo "Start: $(date)"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID}"
echo ""
echo "merged_dge.rds mtime: $(stat -c %y ../results/integration/merged_dge.rds)"
echo ""

# Re-run Script 05 with DREAM_MODEL=base (default) on correct merged_dge.rds
DREAM_MODEL=base Rscript 05_dream_mega_analysis.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 05_dream_mega_analysis.R"; exit 1; fi

echo ""
echo "=== Verifying dream_results.csv ==="
Rscript -e '
library(data.table)
d <- fread("../results/integration/dream_results.csv")
cat("Total genes:", nrow(d), "\n")
cat("DEGs (padj<0.1):", nrow(d[padj < 0.1]), "\n")
cat("DEGs (padj<0.1, |LFC|>0.8):", nrow(d[padj < 0.1 & abs(logFC) > 0.8]), "\n")
' 2>&1

echo ""
echo "=== Script 07: Re-confirm consensus DEGs ==="
# Run Script 07 to ensure consensus_degs.csv matches the restored dream_results.csv
Rscript 07_consensus_degs.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 07_consensus_degs.R"; exit 1; fi

echo ""
echo "=== Restore complete: $(date) ==="
