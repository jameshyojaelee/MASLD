#!/bin/bash
#SBATCH --job-name=nclnasm
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/nclnasm_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/nclnasm_%j.err

set -euo pipefail
umask 027

: "${NONCODING_PROGRAM_LANE_ROOT:?Set NONCODING_PROGRAM_LANE_ROOT}"
: "${NONCODING_PROGRAM_ASSEMBLY:?Set NONCODING_PROGRAM_ASSEMBLY}"

export PYTHONPATH="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/manuscript/noncoding_resource"
/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python \
  scripts/manuscript/noncoding_resource/assemble_program_without_lncrna.py \
  --lane-root "${NONCODING_PROGRAM_LANE_ROOT}" \
  --program-content "RNA-seq/results/noncoding_resource/candidates/noncoding-bulk-2026-08-11/program_annotation/run_19749485/program_lncrna_content.tsv" \
  --output "${NONCODING_PROGRAM_ASSEMBLY}"
