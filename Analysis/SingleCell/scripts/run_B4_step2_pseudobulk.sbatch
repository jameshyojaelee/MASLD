#!/bin/bash
#SBATCH --job-name=B4_step2_pseudobulk
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=16
#SBATCH --mem=512G
#SBATCH --time=48:00:00
#SBATCH --output=Analysis/SingleCell/logs/B4_step2_pseudobulk_%j.out
#SBATCH --error=Analysis/SingleCell/logs/B4_step2_pseudobulk_%j.err

set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p Analysis/SingleCell/logs Analysis/SingleCell/results_gpu_v2/sex_celltype

# Relax `-u` around activate (rapids_singlecell env may have similar binutils hook).
set +u
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rapids_singlecell
set -u

python Analysis/SingleCell/scripts/351_pseudobulk_per_celltype_sex.py \
  --atlas Analysis/SingleCell/results_gpu_v2/integrated_atlas.h5ad \
  --sex Analysis/SingleCell/metadata/snrna_donor_sex.csv \
  --outdir Analysis/SingleCell/results_gpu_v2/sex_celltype \
  --min-cells 10

echo "Step 2 done: $(date)"
