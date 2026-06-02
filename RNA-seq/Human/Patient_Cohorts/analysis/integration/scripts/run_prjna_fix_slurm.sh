#!/bin/bash
# RETIRED 2026-05-15 — PRJNA512027 permanently removed from pipeline.
#
# This sbatch wrapper re-ran 02 -> 03 -> 04 -> 05 -> 07 -> 08 -> 09 -> 10 -> 12
# after applying PRJNA512027-specific normalization fixes (RLE + library_batch
# covariate; dataset_subbatch random effect). The cohort was permanently
# dropped due to the L0/S0 library-prep x disease confound. For the equivalent
# canonical re-run, use the standard integration cascade entry point.
# Kept as a stub for filename provenance.

echo "run_prjna_fix_slurm.sh: RETIRED 2026-05-15 (PRJNA512027 removed). Use run_full_5cohort_slurm.sh instead."
exit 0
