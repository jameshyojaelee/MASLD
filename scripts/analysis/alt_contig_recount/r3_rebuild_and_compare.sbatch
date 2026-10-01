#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=90:00:00
# R3: build orig/ and corr/ run roots (R2), rerun the F_five integration
# (03_integrate_counts.R) and pooled limma-voom QW model (05h) in each with the
# F_five arm settings, then compare. orig/ must reproduce the F_five DEG table.
# Usage: sbatch --job-name=<n> --dependency=afterok:<R1> --output=<out>/logs/%j.out r3_rebuild_and_compare.sbatch <out_dir>
set -euo pipefail
out=$1
repo=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
env=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq
frozen=$repo/RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BULK-F-FIVE/frozen_model_inputs
scripts=$repo/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts
src=$repo/scripts/analysis/alt_contig_recount
gtf=/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz
export PYTHONNOUSERSITE=1

"$env/bin/Rscript" "$src/r2_build_roots.R" "$out/counts" "$gtf" "$out/roots" 2>&1 | tee "$out/r2.log"

for label in orig corr; do
  root=$out/roots/$label
  mkdir -p "$root/provenance"  # 03_integrate_counts.R writes here (line 142) before it creates the dir (line 171)
  MASLD_RUN_ROOT=$root MASLD_INPUT_RUN_ROOT=$root MASLD_COUNTS_MODE=fragment \
  MASLD_QC_MODE=reuse_F_legacy MASLD_FILTER_SCOPE=canonical_five MASLD_GENE_MODE=native \
  MASLD_CONFIG_PATH=$frozen/human_datasets.yaml \
    "$env/bin/Rscript" "$scripts/03_integrate_counts.R" > "$out/$label.03.log" 2>&1
  MASLD_RUN_ROOT=$root MASLD_CONFIG_PATH=$frozen/human_datasets.yaml \
  MASLD_GENE_METADATA_PATH=$frozen/gencode_v49_gene_metadata.tsv.gz \
  MASLD_SANITY_REFERENCE=$repo/results/remediation/bg001/bg001-fragment-v211-gencode49-20260807T195243Z/arms/F_five/results/integration/deg_results.csv \
    "$env/bin/Rscript" "$scripts/05h_limma_voom_qw_canonical.R" > "$out/$label.05h.log" 2>&1
  grep -E "genes=|RAW  ref|CANONICAL TREAT" "$out/$label.05h.log" || true
done

"$env/bin/python" "$src/r4_compare.py" --out "$out" \
  --reference $repo/results/remediation/bg001/bg001-fragment-v211-gencode49-20260807T195243Z/arms/F_five/results/integration/deg_results.csv
