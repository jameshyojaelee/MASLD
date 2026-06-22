#!/bin/bash
#SBATCH --job-name=limma
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/extreme_de_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/extreme_de_%j.err
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=48:00:00

# Extreme-phenotype sensitivity DE (definite-disease vs strict-control, 05i)
# + the default (n=200) paired heatmap. cpu partition, default nslab QOS.
set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
set +u; micromamba activate rnaseq; set -u

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$ROOT"
echo "=== extreme-phenotype DE started: $(date) | job $SLURM_JOB_ID ==="

echo "--- 05i: definite-vs-strict-control LVQW C2 ---"
Rscript RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/05i_extreme_phenotype_de.R 2>&1

echo "--- default paired heatmap (n=200, no suffix) ---"
TOP_N_DEG=200 HEATMAP_OUT_SUFFIX="" \
  Rscript scripts/figures/figS_extreme_phenotype_heatmap.R 2>&1

echo "=== done: $(date) ==="
