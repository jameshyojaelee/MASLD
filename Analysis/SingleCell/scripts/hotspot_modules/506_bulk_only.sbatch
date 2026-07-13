#!/bin/bash
# Re-run only 506 bulk-replication (R/rnaseq) on the current post-remediation
# module_genes.tsv + the canonical C2 limma-voom-qw DEG. Independent of the LOO
# array; can run immediately. (507 olink/spatial is skipped — inputs are absent.)
#SBATCH --job-name=hotspot_506
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=Analysis/SingleCell/scripts/hotspot_modules/logs/506_bulk_%j.out
#SBATCH --error=Analysis/SingleCell/scripts/hotspot_modules/logs/506_bulk_%j.err
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p Analysis/SingleCell/scripts/hotspot_modules/logs
micromamba run -n rnaseq Rscript Analysis/SingleCell/scripts/hotspot_modules/506_bulk_replication.R
