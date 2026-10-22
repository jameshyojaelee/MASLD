#!/bin/bash
#SBATCH --job-name=pazoki_orient_check
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/orientation_checks/logs/slurm-%j.out

# Step 1: is the effect/other allele labelling of the Pazoki 2022 PDFF file (GCST90267352) reversed?
# Read-only comparison of the as-downloaded file and its three derived copies against two
# independently oriented PDFF GWAS. Writes a new timestamped directory; modifies no input.
set -euo pipefail
REPO=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
SC=$REPO/scripts/analysis/gwas_orientation
PY=/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python
TS=$(date -u +%Y%m%dT%H%M%SZ)
OUT=$REPO/GWAS/orientation_checks/pazoki_pdff-$TS
mkdir "$OUT"   # fails if it exists: never overwrite an earlier run
echo "host $(hostname)  start $(date -u +%FT%TZ)  job ${SLURM_JOB_ID:-NA}  out $OUT"
sha256sum "$SC/pazoki_orientation_check.py" "$SC/run_orientation_check.sbatch" | tee "$OUT/script_sha256.txt"
export PYTHONHASHSEED=0 OMP_NUM_THREADS=4
"$PY" -I "$SC/pazoki_orientation_check.py" "$OUT"
echo "end $(date -u +%FT%TZ)"
