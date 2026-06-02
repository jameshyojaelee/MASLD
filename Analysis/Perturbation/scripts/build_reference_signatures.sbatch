#!/bin/bash
#SBATCH --job-name=d2_ref_curator
#SBATCH --partition=gpu
#SBATCH --gres=gpu:0
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=Analysis/Perturbation/scripts/logs/d2_ref_curator_%j.out
#SBATCH --error=Analysis/Perturbation/scripts/logs/d2_ref_curator_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

mkdir -p Analysis/Perturbation/scripts/logs

micromamba run -n rnaseq python Analysis/Perturbation/scripts/build_reference_signatures.py
