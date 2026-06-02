#!/bin/bash
#SBATCH --job-name=ora_fix_v3
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=logs/ora_fix_v3_%j.log

set -euo pipefail
BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${BASE}"
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
export MASLD_PROJECT_ROOT="${BASE}"
echo "=== ORA fix v3: $(date) ==="
Rscript Analysis/Cross_Species_Concordance/scripts/02c_ora_concordance.R
echo "=== Completed: $(date) ==="
