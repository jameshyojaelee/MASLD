#!/bin/bash
#SBATCH --job-name=ggplot2
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=scripts/figures/logs/convergence_rerender_%j.out
#SBATCH --error=scripts/figures/logs/convergence_rerender_%j.err

set -o pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
export MASLD_PROJECT_ROOT="${BASE}"
mkdir -p "${BASE}/scripts/figures/logs"

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

cd "${BASE}/scripts/figures"

echo "[$(date '+%H:%M:%S')] figS_convergence.R"
Rscript figS_convergence.R || echo "WARN figS_convergence.R failed"

echo "[$(date '+%H:%M:%S')] figS_evidence_convergence.R"
Rscript figS_evidence_convergence.R || echo "WARN figS_evidence_convergence.R failed"

echo "[$(date '+%H:%M:%S')] figS_convergence_evidence.R"
Rscript figS_convergence_evidence.R || echo "WARN figS_convergence_evidence.R failed"

echo "[$(date '+%H:%M:%S')] DONE"
