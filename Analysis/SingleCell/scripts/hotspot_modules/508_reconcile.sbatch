#!/bin/bash
# Backup the canonical multi-evidence atlas, then re-run 508 to rebuild
# all_modules.tsv + the 6 hotspot_* atlas columns from the now-coherent
# post-remediation inputs (504 novelty, 503b LOO, 506 bulk, 505 phenotype,
# current module_genes). Chained afterok on 503b + 504 + 506.
#SBATCH --job-name=hotspot_508
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=Analysis/SingleCell/scripts/hotspot_modules/logs/508_reconcile_%j.out
#SBATCH --error=Analysis/SingleCell/scripts/hotspot_modules/logs/508_reconcile_%j.err
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

ATLAS=RNA-seq/results/multi_evidence/multi_evidence_atlas.csv
BAK="${ATLAS}.bak_pre_hotspot_reconcile_20260630"
if [ ! -f "$BAK" ]; then
  cp -p "$ATLAS" "$BAK"
  echo "[508_reconcile] backed up atlas -> $BAK"
else
  echo "[508_reconcile] backup already exists: $BAK (not overwriting)"
fi

micromamba run -n rnaseq Rscript Analysis/SingleCell/scripts/hotspot_modules/508_atlas_integration.R
