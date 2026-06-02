#!/bin/bash
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -uo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
for s in figS_susie_gwas_volcano.R figS_susie_expression_scatters.R figS_susie_celltype_coloc_heatmap.R; do
  echo "=== Running $s ==="
  Rscript scripts/figures/$s 2>&1 || echo "FAIL: $s"
done
echo "=== Done ==="
