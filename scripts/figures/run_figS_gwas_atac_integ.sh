#!/bin/bash
# micromamba activate must come BEFORE set -euo pipefail (ADDR2LINE bug in activate-binutils)
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
for s in figS05_gwas_atac.R fig_gwas_atac_regulons.R figS_twas_coloc_intact.R figS_gwas_celltype.R; do
  echo "=== Running $s ==="
  Rscript scripts/figures/$s 2>&1 || echo "FAIL: $s"
done
