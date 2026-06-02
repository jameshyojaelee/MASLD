#!/bin/bash
#SBATCH --partition=cpu --cpus-per-task=4 --mem=64G --time=48:00:00
#SBATCH --job-name=masld_kcommit --output=RNA-seq/logs/masld_kcommit_%j.out --error=RNA-seq/logs/masld_kcommit_%j.err
set -eo pipefail
mkdir -p RNA-seq/logs
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

# Phase 2: Rewrite Script 44 output at chosen k
echo "=== Phase 2: Script 44 (k-program subtyping) ==="
Rscript RNA-seq/44_molecular_subtyping.R

# Phase 3a: Script 208 multinomial COLOC
echo ""
echo "=== Phase 3a: Script 208 (k-program COLOC) ==="
Rscript RNA-seq/208_subtype_coloc.R

# Phase 3b: Script 217 atlas merge
echo ""
echo "=== Phase 3b: Script 217 (atlas columns) ==="
Rscript RNA-seq/217_stratified_causal_atlas.R

# Phase 4: Sanity checks
echo ""
echo "=== Phase 4: Sanity ==="
Rscript RNA-seq/93_kprogram_sanity.R

# Phase 3c: figS08 regeneration
echo ""
echo "=== Phase 3c: figS08 ==="
Rscript scripts/figures/figS08_nmf_subtyping.R

echo ""
echo "=== k-program commit pipeline complete ==="
