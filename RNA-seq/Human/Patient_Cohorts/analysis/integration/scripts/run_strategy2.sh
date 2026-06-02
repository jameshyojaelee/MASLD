#!/bin/bash
#SBATCH --job-name=liver_deconv_attrib
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/causal_inference/logs/deconv_attrib_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/causal_inference/logs/deconv_attrib_%j.err
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=32
#SBATCH --mem=300G
#SBATCH --time=06:00:00

set -euo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration

echo "=== Strategy 2: Deconvolution Attribution Analysis ==="
echo "Started: $(date)"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo ""

Rscript scripts/25_deconv_attribution.R 2>&1
if [ $? -ne 0 ]; then echo "FAILED: 25_deconv_attribution.R"; exit 1; fi

echo ""
echo "=== Strategy 2 complete: $(date) ==="
