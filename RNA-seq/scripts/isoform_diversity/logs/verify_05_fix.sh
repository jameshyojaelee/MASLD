#!/bin/bash
#SBATCH --job-name=integrate
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/verify_05_fix_%j.out
#SBATCH --error=logs/verify_05_fix_%j.err
# Verify the 05_integrate.R stratum fix: (1) byte-stability across two runs (determinism),
# (2) quadrant/dtu_pos diff OLD(pooled, .prepatch_pooled_master.tsv) vs NEW(total-RNA).
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/scripts/isoform_diversity
RD=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/isoform_diversity/human

echo "=========== RUN 1 (new default stratum = group_total-RNA) ==========="
micromamba run -n dtu Rscript 05_integrate.R human
cp "$RD/isoform_master_human.tsv" "$RD/.verify_run1.tsv"

echo "=========== RUN 2 (determinism) ==========="
micromamba run -n dtu Rscript 05_integrate.R human
cp "$RD/isoform_master_human.tsv" "$RD/.verify_run2.tsv"

echo "=========== MD5 (run1 vs run2 must match) ==========="
md5sum "$RD/.verify_run1.tsv" "$RD/.verify_run2.tsv"

echo "=========== quadrant OLD(pooled) vs NEW(total-RNA) ==========="
micromamba run -n dtu Rscript -e '
library(data.table)
RD <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/isoform_diversity/human"
o <- fread(file.path(RD, ".prepatch_pooled_master.tsv"))
n <- fread(file.path(RD, "isoform_master_human.tsv"))
cat("dims OLD/NEW:", nrow(o),"x",ncol(o)," | ", nrow(n),"x",ncol(n), "\n\n")
cat("OLD quadrant (pooled dtu_human_group):\n"); print(table(o$quadrant, useNA="ifany"))
cat("  dtu_pos:", sum(o$dtu_pos,na.rm=TRUE), " pure_switch_strong:", sum(o$pure_switch_strong,na.rm=TRUE), "\n\n")
cat("NEW quadrant (total-RNA):\n"); print(table(n$quadrant, useNA="ifany"))
cat("  dtu_pos:", sum(n$dtu_pos,na.rm=TRUE), " pure_switch_strong:", sum(n$pure_switch_strong,na.rm=TRUE), "\n\n")
cat("disease_iso non-NA OLD/NEW:", sum(!is.na(o$disease_iso)), "/", sum(!is.na(n$disease_iso)), "\n")
'
echo "=========== cleanup temp run copies ==========="
rm -f "$RD/.verify_run1.tsv" "$RD/.verify_run2.tsv"
echo "=========== DONE ==========="
