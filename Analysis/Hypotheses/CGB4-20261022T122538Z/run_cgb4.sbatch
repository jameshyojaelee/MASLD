#!/bin/bash
#SBATCH --job-name=cgb4
#SBATCH --account=nslab
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=3
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CGB4-20261022T122538Z/logs/slurm-%x-%j.out
set -euo pipefail
MRUN=$(ls -d /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CG0-20261022T122300Z/magma/run-* | sort | tail -1)
OUT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CGB4-20261022T122538Z/out-$(date -u +%Y%m%dT%H%M%SZ)
bash /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CGB4-20261022T122538Z/scripts/cgb4_magma_properties.sh /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CG0-20261022T122300Z/out-20261022T122649Z $MRUN $OUT
echo "$MRUN" > $OUT/magma_run_used.txt
sha256sum /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CGB4-20261022T122538Z/scripts/* /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CGB4-20261022T122538Z/PRESPEC_CGB4.md > $OUT/script_sha256.txt
