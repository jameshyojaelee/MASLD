#!/bin/bash
# Generic R script runner for SLURM.
# Usage: sbatch [slurm opts] run_rscript.sh script1.R [script2.R ...]
# Runs each R script sequentially within one job.

set -euo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts

for script in "$@"; do
    name=$(basename "$script" .R)
    echo "=== ${name} ==="
    echo "Started: $(date)"
    Rscript "${script}" 2>&1
    echo "Completed: $(date)"
    echo ""
done
