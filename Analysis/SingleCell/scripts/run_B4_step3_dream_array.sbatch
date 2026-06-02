#!/bin/bash
#SBATCH --job-name=B4_step3_dream
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --array=0-4
#SBATCH --output=Analysis/SingleCell/logs/B4_step3_dream_%A_%a.out
#SBATCH --error=Analysis/SingleCell/logs/B4_step3_dream_%A_%a.err

set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p Analysis/SingleCell/logs

# rnaseq env's binutils activate hook references unbound ADDR2LINE; relax `-u` around activate.
set +u
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

CTS=(hepatocytes macrophages fibroblasts endothelial cholangiocytes)
CT="${CTS[$SLURM_ARRAY_TASK_ID]}"
echo "Cell type: $CT"

Rscript Analysis/SingleCell/scripts/352_dream_per_celltype_sex.R "$CT"

echo "Step 3 ($CT) done: $(date)"
