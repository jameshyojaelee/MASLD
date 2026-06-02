#!/bin/bash
#SBATCH --job-name=deconv
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/scripts/featurecounts_s2/logs/deconv_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/scripts/featurecounts_s2/logs/deconv_%j.err

set -eo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

set -u

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

echo "=== 25: deconvolution attribution (-s 2; reuses existing MuSiC scRNA-ref props). Writes results/causal_inference/deconv_attribution_scores.csv ==="
echo "Started: $(date)  Host: $(hostname)  CPUs: ${SLURM_CPUS_PER_TASK}"

Rscript RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/25_deconv_attribution.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 25_deconv_attribution.R"; exit 1; fi

echo "Finished: $(date)"
