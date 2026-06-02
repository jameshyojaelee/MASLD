#!/bin/bash
#SBATCH --job-name=liver_dream_variants
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=32
#SBATCH --mem=400G
#SBATCH --time=12:00:00
#SBATCH --output=logs/dream_variants_%j.out
#SBATCH --error=logs/dream_variants_%j.err

SCRIPTS_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts"
cd "$SCRIPTS_DIR"
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

echo "=== Dream Model Variant Comparison ==="
echo "Start: $(date)"

for MODEL in base randslope robust randslope_robust; do
  echo ""
  echo "--- Running model: $MODEL ---"
  DREAM_MODEL=$MODEL Rscript 05_dream_mega_analysis.R
done

# Also run 5-cohort sensitivity (exclude GSE213621)
echo ""
echo "--- Running 5-cohort sensitivity (exclude GSE213621) ---"
DREAM_MODEL=base EXCL_EXTRA=GSE213621 Rscript 05_dream_mega_analysis.R

echo ""
echo "=== Comparing results ==="
Rscript -e '
library(data.table)
res_dir <- "../results/integration"
for (f in list.files(res_dir, pattern="dream_results.*\\.csv$", full.names=TRUE)) {
  d <- fread(f)
  sig <- d[padj < 0.1 & abs(logFC) > 0.8]
  cat(basename(f), ": ", nrow(sig), " DEGs (padj<0.1, |LFC|>0.8),",
      " up=", sum(sig$logFC > 0), " down=", sum(sig$logFC < 0), "\n")
}
'

echo "Done: $(date)"
