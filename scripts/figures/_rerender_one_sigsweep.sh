#!/bin/bash
#SBATCH --job-name=figregen
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/figures/logs/rerender_sigsweep_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/figures/logs/rerender_sigsweep_%j.err
set -o pipefail
BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
export MASLD_PROJECT_ROOT="$BASE"; cd "$BASE"
eval "$(micromamba shell hook --shell bash)"; micromamba activate rnaseq
echo ">>> figS_disease_signature_sweep.R (no per-script timeout; rely on 48h wall)"
t0=$SECONDS
Rscript scripts/figures/figS_disease_signature_sweep.R
rc=$?
echo "[exit=$rc] elapsed $((SECONDS-t0))s"
