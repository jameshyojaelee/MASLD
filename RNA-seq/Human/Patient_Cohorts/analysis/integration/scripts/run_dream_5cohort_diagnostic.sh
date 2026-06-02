#!/bin/bash
#SBATCH --job-name=liver_dream_diag
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=32
#SBATCH --mem=300G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/dream_5cohort_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/dream_5cohort_%j.err
#SBATCH --export=ALL

export EXCL_EXTRA=GSE213621

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
set +u; micromamba activate rnaseq; set -u
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts

echo "5-cohort dream diagnostic (excl GSE167523 + GSE213621): $(date)"
Rscript analysis/integration/scripts/05_dream_mega_analysis.R 2>&1
echo "Script 05 complete: $(date)"
Rscript analysis/integration/scripts/07_consensus_degs.R 2>&1
echo "Script 07 complete: $(date)"
echo "Done: $(date)"
