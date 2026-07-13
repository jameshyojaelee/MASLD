#!/bin/bash
#SBATCH --job-name=decima
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/decima_celltype_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/decima_celltype_%j.log

set -euo pipefail
export MAMBA_ROOT_PREFIX=/gpfs/commons/home/jameslee/micromamba
export PYTHONNOUSERSITE=1                       # isolate from ~/.local leaked deps
SF=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc
export HF_HOME=$SF/hf_cache

DPY=/gpfs/commons/home/jameslee/micromamba/envs/decima/bin/python
echo "host $(hostname) | $(date)"
nvidia-smi -L
$DPY -c "import torch; print('cuda avail', torch.cuda.is_available(), torch.cuda.get_device_name(0))"
$DPY /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/src/76_decima_celltype.py
echo "EXIT $? | $(date)"
