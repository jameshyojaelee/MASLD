#!/bin/bash
#SBATCH --job-name=cas13_seq_noov
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=12
#SBATCH --mem=120G
#SBATCH --time=08:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Cas13_Library_Design/scripts/guides/guide_overlap_fix_20260713/slurm-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Cas13_Library_Design/scripts/guides/guide_overlap_fix_20260713/slurm-%j.err
set -eo pipefail

PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
RUN=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Cas13_Library_Design/scripts/guides/guide_overlap_fix_20260713

cd "$RUN"
"$PY" validate_submission.py
"$PY" run_isolated_build.py
"$PY" combine_isolated.py
"$PY" audit_library_overlap.py output/library_sequence_nonoverlap.csv \
  --report output/rebuilt_overlap_audit.csv
