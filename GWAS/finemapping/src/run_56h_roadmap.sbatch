#!/bin/bash
#SBATCH --job-name=roadmap_e066
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/56h_roadmap_%j.out
#SBATCH --error=logs/56h_roadmap_%j.err

set -eo pipefail

PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$PROJ"
mkdir -p data/external/roadmap_liver_E066 \
         GWAS/finemapping/src/logs \
         GWAS/finemapping/results/gwas_atac

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

echo "[$(date)] Starting Roadmap E066 replication pipeline (job $SLURM_JOB_ID)"
echo "Host: $(hostname)"
echo "CPUs: $SLURM_CPUS_PER_TASK"

# Step 1: download + liftover (idempotent)
if [[ ! -s data/external/roadmap_liver_E066/E066_DNase_hg38.broadPeak.gz ]]; then
    bash GWAS/finemapping/src/56h_download_roadmap_E066.sh
else
    echo "Roadmap files already present, skipping download"
fi

# Step 2: replication analysis
Rscript GWAS/finemapping/src/56h_roadmap_replication.R

echo "[$(date)] Done."
