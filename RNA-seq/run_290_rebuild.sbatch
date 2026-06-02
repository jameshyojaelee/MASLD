#!/bin/bash
#SBATCH --job-name=290_atlas_rebuild
#SBATCH --partition=cpu
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --output=logs/290_atlas_rebuild_%j.out
#SBATCH --error=logs/290_atlas_rebuild_%j.err

set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
micromamba run -n spatial python 290_edge_annotation_atlas.py
