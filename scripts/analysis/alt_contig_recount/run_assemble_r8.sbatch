#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=90:00:00
#SBATCH --job-name=assemble_r8
#SBATCH --output=logs/alt_contig_recount/assemble_r8-%j.out
# Assemble both count arms for the 774-unit and 738-unit (strict depth) rosters.
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
candidate=$PWD/figures/candidates/resource-corrections-20261005-v1
counts=$candidate/verified_counts_r8_20261007
n=$(ls "$counts"/samples/*/count_completion.json | wc -l)
[[ "$n" == 776 ]] || { echo "Only $n of 776 libraries counted" >&2; exit 1; }
rscript=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
for roster in 774:identity_resource_qualified 738:identity_resource_strict_depth; do
  for arm in primary single_candidate_credit; do
    "$rscript" --vanilla scripts/analysis/alt_contig_recount/assemble_verified_counts.R \
      "$candidate/annotation_resolution_r8/resolution_summary.json" \
      "$candidate/${roster#*:}/bulk_primary_sample_ids.tsv" "$counts/samples" \
      "$counts/assembled_${roster%%:*}_${arm}" "$arm"
  done
done
