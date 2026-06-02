#!/bin/bash
#SBATCH --job-name=160b_fix_hcc
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=results/prognosis_v2/160b_fix_hcc_%j.out
#SBATCH --error=results/prognosis_v2/160b_fix_hcc_%j.err

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration

echo "=== 160b: Fix HCC molecular score gene mapping ==="
echo "Started: $(date)"
echo "Node: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"

micromamba run -n spatial python scripts/160b_fix_hcc_score.py

echo "Finished: $(date)"
