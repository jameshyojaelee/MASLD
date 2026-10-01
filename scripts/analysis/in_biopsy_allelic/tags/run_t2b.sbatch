#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=24:00:00
# Spec v1 tag table (1.2), then B-T2 per-cohort het rule, phi and the Gate A aggregate inputs.
# Usage: sbatch --job-name=<n> --output=<out>/logs/%j.out run_t2b.sbatch <out>
set -euo pipefail
REPO=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
EX=$REPO/Analysis/MASLD_Model_Benchmark/executions
RS=$REPO/RNA-seq/Human/Patient_Cohorts/results/rna_genotypes_restricted
PY="env PYTHONNOUSERSITE=1 /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python"
T=$REPO/scripts/analysis/in_biopsy_allelic/tags
QC=$EX/in-biopsy-allelic-obs-20260926T155651Z/qc/v2-20260929T143503Z
AT=$EX/in-biopsy-allelic-tags-v2-20260929T134419Z/anc_tags-20260929T200201Z
$PY $T/build_tag_table.py \
  --t3prime $EX/in-biopsy-allelic-obs-20260926T155651Z/t3prime/v2-20260929T200244Z/tag_exclusions_v2.tsv \
  --freq $AT/tag_2pq_by_superpop.tsv --ld $AT/tag_ld_by_superpop.tsv \
  --ffpe-drop GSE213621,GSE240729 --out "$1/tag_table"
$PY $T/t2b_cohort_rule.py \
  --geuvadis-work $EX/in-biopsy-allelic-tags-20260925T194248Z/t2_work \
  --crosswalk $EX/in-biopsy-allelic-a1-full-20260924T025858Z/seal-20260929T141204Z/frozen_crosswalk.tsv \
  --ad $RS/tags-ad-v2-20260929T142228Z/tag_allele_depths.tsv.gz \
  --tag-table "$1/tag_table/tag_table.tsv" --null-bias $QC/cohort_null_site_ref_bias.tsv \
  --error-rate $QC/cohort_error_rate.tsv --tag-freq $AT/tag_2pq_by_superpop.tsv \
  --e-calibrated $EX/in-biopsy-allelic-obs-20260926T155651Z/qc/r1c-20260930T123747Z/cohort_model_inputs.tsv \
  --ancestry $RS/obs-20260926T155651Z/ancestry/full-20260929T142807Z/library_ancestry_pcs.tsv \
  --out "$1"
$PY -m pip freeze > "$1/pip_freeze.txt"
sha256sum $T/build_tag_table.py $T/t2b_cohort_rule.py > "$1/sha256_scripts.txt"
