#!/usr/bin/env bash
#SBATCH --job-name=anndata
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/build_hmsma_anndata_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/logs/build_hmsma_anndata_%j.out
set -o pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
micromamba run -n spatial python \
  Analysis/Spatial/scripts/53_build_hmsma_anndata.py --min-umi 500
echo "EXIT=$?"
