#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --job-name=343j_dubHealthy
#SBATCH --output=Analysis/SingleCell/logs/343j_dubHealthy_%j.out
#SBATCH --error=Analysis/SingleCell/logs/343j_dubHealthy_%j.err
set -eo pipefail
mkdir -p Analysis/SingleCell/logs
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
micromamba run -n rapids_singlecell python \
    Analysis/SingleCell/scripts/343j_dubious_healthy_audit.py
