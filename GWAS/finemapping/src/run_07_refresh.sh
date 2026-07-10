#!/bin/bash
#SBATCH --job-name=07_refresh
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/src/logs/07_refresh_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/src/logs/07_refresh_%j.err
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping
# 07 -> 07b -> 07 : the first pass writes susie_coloc_all_gwas.csv; 07b rebuilds
# ld_contamination_clusters.csv over the full (MVP-inclusive) master; the second
# pass merges the fresh LD-contamination flags into gene_level_coloc.csv.
# (Without the 07b step, 07 merges a STALE cluster file and any new MVP
# LD-contaminated locus gets ld_cluster_flag=NA -> miscounted as independent
# gene-level evidence.)
Rscript src/07_combine_susie_coloc.R
Rscript src/07b_flag_ld_clusters.R
Rscript src/07_combine_susie_coloc.R
