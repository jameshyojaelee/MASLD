#!/bin/bash
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
for s in figS04_coloc.R figS09_multi_ancestry_coloc.R figS09_locus_zoom.R figS09d_mesusie_susiex_comparison.R figS_finemapping_concordance.R; do
  echo "=== Running $s ==="
  Rscript scripts/figures/$s 2>&1 || echo "FAIL: $s"
done
