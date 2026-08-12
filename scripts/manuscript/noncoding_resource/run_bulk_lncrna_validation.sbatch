#!/bin/bash
#SBATCH --job-name=nclnval
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/nclnval_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/nclnval_%j.err

set -euo pipefail
umask 027

: "${NONCODING_BULK_CANDIDATE:?Set NONCODING_BULK_CANDIDATE}"
: "${NONCODING_BULK_REVIEW:?Set NONCODING_BULK_REVIEW}"

python3 scripts/manuscript/noncoding_resource/validate_bulk_lncrna.py \
  --candidate "${NONCODING_BULK_CANDIDATE}" \
  --output "${NONCODING_BULK_REVIEW}"
