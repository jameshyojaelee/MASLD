#!/bin/bash
#SBATCH --job-name=stg118_bayesprism
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=210G
#SBATCH --time=72:00:00
#SBATCH --output=scripts/logs/118_bayesprism_%j.out
#SBATCH --error=scripts/logs/118_bayesprism_%j.err

set -euo pipefail
source ~/.bashrc
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration

# Part A: Prepare scRNA reference (Python, rapids_singlecell env)
REF_DIR="results/progression/scrna_reference"
if [ -f "${REF_DIR}/scrna_reference_summary.csv" ] && [ -f "${REF_DIR}/scrna_counts.csv.gz" ]; then
    echo "=== Part A: scRNA reference already exists, skipping ==="
else
    echo "=== Part A: Preparing scRNA reference ==="
    micromamba activate rapids_singlecell
    python scripts/118a_prepare_scrna_reference.py
    micromamba deactivate
fi

# Part B: Run BayesPrism deconvolution (R, rnaseq env)
echo "=== Part B: Running BayesPrism ==="
micromamba activate rnaseq
Rscript scripts/118b_bayesprism_hires.R

echo "=== 118: ALL COMPLETE ==="
