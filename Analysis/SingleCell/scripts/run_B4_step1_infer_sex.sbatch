#!/bin/bash
#SBATCH --job-name=B4_step1_inferSex
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=Analysis/SingleCell/logs/B4_step1_inferSex_%j.out
#SBATCH --error=Analysis/SingleCell/logs/B4_step1_inferSex_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p Analysis/SingleCell/logs Analysis/SingleCell/metadata

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rapids_singlecell

python Analysis/SingleCell/scripts/350_infer_sex_snrna.py \
  --atlas Analysis/SingleCell/results_gpu_v2/integrated_atlas.h5ad \
  --multiome-ref Analysis/ATAC/Human_Multiome/results/sex_validation/donor_sex_inference.csv \
  --out Analysis/SingleCell/metadata/snrna_donor_sex.csv \
  --summary Analysis/SingleCell/metadata/snrna_donor_sex_summary.json

echo "Step 1 done: $(date)"
