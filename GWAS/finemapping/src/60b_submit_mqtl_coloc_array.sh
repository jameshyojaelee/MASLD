#!/bin/bash
# 60b_submit_mqtl_coloc_array.sh
# Submit the metabolite/lipid-species COLOC array (60b_mqtl_coloc.R).
# Designed to run AFTER 60a_download_mqtl_sumstats.sh finishes.
# Usage: bash 60b_submit_mqtl_coloc_array.sh [--dependency=afterok:<jobid>]
#
# Two arms submitted independently:
#   Chen 2023  — 142 metabolite files -> N_SHARDS=142 (1 shard per file, ABF-only)
#   Ottensmann — N_SHARDS based on ottensmann accession map

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SRC="$BASE/GWAS/finemapping/src"
LOGS="$BASE/GWAS/finemapping/logs"
mkdir -p "$LOGS"

DEPENDENCY="${1:-}"          # e.g. --dependency=afterok:16963028

# ---- count shards per arm ----
CHEN_LIST="$BASE/data/external/chen2023_mqtl/chen2023_biomarker_download_list.tsv"
LIPID_LIST="$BASE/data/external/lipidqtl/ottensmann2023_accession_trait_map.tsv"
N_CHEN=$(( $(wc -l < "$CHEN_LIST") - 1 ))   # subtract header
N_LIPID=$(( $(wc -l < "$LIPID_LIST") - 1 ))

echo "Chen metabolites:  $N_CHEN"
echo "Ottensmann lipids: $N_LIPID"
echo "Dependency:        ${DEPENDENCY:-none}"

RBIN="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"

# ---- Chen arm array ----
JOB_CHEN=$(sbatch --parsable \
  --job-name=coloc \
  --partition=cpu \
  --qos=nslab \
  --array=0-$((N_CHEN-1)) \
  --cpus-per-task=2 \
  --mem=16G \
  --time=48:00:00 \
  --output="$LOGS/60b_chen_%A_%a.out" \
  --error="$LOGS/60b_chen_%A_%a.err" \
  ${DEPENDENCY:+"$DEPENDENCY"} \
  --wrap="export QTL_SOURCE=chen SHARD_ID=\$SLURM_ARRAY_TASK_ID N_SHARDS=$N_CHEN && \
          $RBIN $SRC/60b_mqtl_coloc.R")
echo "Chen COLOC array submitted: $JOB_CHEN"

# ---- Ottensmann arm array ----
JOB_LIPID=$(sbatch --parsable \
  --job-name=coloc \
  --partition=cpu \
  --qos=nslab \
  --array=0-$((N_LIPID-1)) \
  --cpus-per-task=2 \
  --mem=16G \
  --time=48:00:00 \
  --output="$LOGS/60b_lipid_%A_%a.out" \
  --error="$LOGS/60b_lipid_%A_%a.err" \
  ${DEPENDENCY:+"$DEPENDENCY"} \
  --wrap="export QTL_SOURCE=ottensmann SHARD_ID=\$SLURM_ARRAY_TASK_ID N_SHARDS=$N_LIPID && \
          $RBIN $SRC/60b_mqtl_coloc.R")
echo "Lipid  COLOC array submitted: $JOB_LIPID"

echo ""
echo "Monitor: squeue --name=coloc -u \$USER"
echo "Combine after: $RBIN $SRC/61b_combine_metabolite_coloc.R"
