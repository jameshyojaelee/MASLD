#!/bin/bash
# Submit the stratified diversity + DTU fan-out for one species.
# Assumes 02_aggregate (tx matrices + samples.tsv) is already complete.
# Pooled diversity is chained separately; this adds per-chemistry diversity + all DTU.
#   bash run_fanout.sh <human|mouse> [dtu]   # pass "dtu" to also submit DTU jobs
set -euo pipefail
SP=$1; WITH_DTU=${2:-}
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/scripts/isoform_diversity

echo "[fanout] $SP diversity (per-chemistry)"
for LT in polyA total-RNA; do
  sbatch --parsable --export=ALL,SPECIES=$SP,LIBTYPE=$LT 03_diversity.sh | sed "s/^/  div $SP $LT job /"
done

if [ "$WITH_DTU" = "dtu" ]; then
  echo "[fanout] $SP DTU group (pooled + per-chemistry)"
  sbatch --parsable --export=ALL,SPECIES=$SP,CONTRAST=group 04_dtu.sh | sed "s/^/  dtu $SP group(pooled) job /"
  for LT in polyA total-RNA; do
    sbatch --parsable --export=ALL,SPECIES=$SP,CONTRAST=group,LIBTYPE=$LT 04_dtu.sh | sed "s/^/  dtu $SP group($LT) job /"
  done
  if [ "$SP" = human ]; then
    echo "[fanout] human DTU fibrosis-stage axis"
    sbatch --parsable --export=ALL,SPECIES=human,CONTRAST=stage 04_dtu.sh | sed "s/^/  dtu human stage job /"
  fi
fi
echo "[fanout] $SP submitted."
