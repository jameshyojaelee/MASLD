#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=6:00:00
# Gate A design export (depth, stage, z; no allele fraction). Usage: sbatch ... run_gateA_export.sbatch <out>
set -euo pipefail
REPO=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
EX=$REPO/Analysis/MASLD_Model_Benchmark/executions
RS=$REPO/RNA-seq/Human/Patient_Cohorts/results/rna_genotypes_restricted
T2B=$EX/in-biopsy-allelic-t2b-20260929T220019Z
PY="env PYTHONNOUSERSITE=1 /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python"
$PY $REPO/scripts/analysis/in_biopsy_allelic/tags/gateA_design_export.py \
  --crosswalk $EX/in-biopsy-allelic-a1-full-20260924T025858Z/seal-20260929T141204Z/frozen_crosswalk.tsv \
  --ad $RS/tags-ad-v2-20260929T142228Z/tag_allele_depths.tsv.gz \
  --tag-table $T2B/tag_table/tag_table.tsv \
  --tag-freq $EX/in-biopsy-allelic-tags-v2-20260929T134419Z/anc_tags-20260929T200201Z/tag_2pq_by_superpop.tsv \
  --ancestry $RS/obs-20260926T155651Z/ancestry/full-20260929T142807Z/library_ancestry_pcs.tsv \
  --library-qc $RS/obs-20260926T155651Z/qc/r1-20260930T121200Z/library_model_inputs.tsv \
  --dup $RS/obs-20260926T155651Z/dup/gene-full-20260929T142521Z/collect/dup_per_library_v2.tsv \
  --metadata $REPO/RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv \
  --kinship $EX/in-biopsy-allelic-a1-full-20260924T025858Z/identity/kinship_related_pairs.tsv \
  --rules $T2B/cohort_rules.tsv --out "$1"
sha256sum $REPO/scripts/analysis/in_biopsy_allelic/tags/gateA_design_export.py > "$1/sha256_script.txt"
