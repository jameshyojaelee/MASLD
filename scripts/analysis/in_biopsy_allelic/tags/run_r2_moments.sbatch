#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=4:00:00
# Conditional-lambda covariate moments and fixed P2 resamples of F. Usage: sbatch ... run_r2_moments.sbatch <out>
set -euo pipefail
REPO=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
EX=$REPO/Analysis/MASLD_Model_Benchmark/executions
PY="env PYTHONNOUSERSITE=1 /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python"
$PY $REPO/scripts/analysis/in_biopsy_allelic/tags/gateA_r2_moments.py \
  --individuals $EX/model-a-gateA-design-20260930T121812Z/design_individuals.tsv \
  --crosswalk $EX/in-biopsy-allelic-a1-full-20260924T025858Z/seal-20260929T141204Z/frozen_crosswalk.tsv \
  --ancestry $REPO/RNA-seq/Human/Patient_Cohorts/results/rna_genotypes_restricted/obs-20260926T155651Z/ancestry/full-20260929T142807Z/library_ancestry_pcs.tsv \
  --draws 2000 --out "$1"
sha256sum $REPO/scripts/analysis/in_biopsy_allelic/tags/gateA_r2_moments.py > "$1/sha256_script.txt"
