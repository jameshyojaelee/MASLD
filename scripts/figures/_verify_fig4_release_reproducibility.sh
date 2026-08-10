#!/usr/bin/env bash
#SBATCH --job-name=data.table
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/fig4_release_repro_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/fig4_release_repro_%j.out

# READ-ONLY verification: can the frozen 2026-07-15-r2 evidence-class validation
# endpoints still be reproduced from the LIVE multi_evidence_atlas.csv?
#
# The release builder (scripts/manuscript/build_evidence_class_release.R:262-281)
# takes its three endpoints from atlas columns best_protein_padj / spatial_is_svg /
# sc_best_padj, joined to the release's own evidence_class_table.tsv on symbol.
# This reproduces that join and prints the resulting denominators/positives so they
# can be compared against evidence_class_validation_summary.tsv.
#
# Writes NOTHING. Does not touch the frozen release.

set -eo pipefail
RBIN=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
echo "host=$(hostname) start=$(date)"

"$RBIN" -e '
suppressPackageStartupMessages(library(data.table))
setDTthreads(4)
B <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
REL <- file.path(B,"RNA-seq/results/manuscript_release/2026-07-15-r2")

atlas <- fread(file.path(B,"RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"),
  select=c("human_symbol","best_protein_padj","spatial_is_svg","sc_best_padj"),
  showProgress=FALSE)
setnames(atlas,"human_symbol","symbol")
atlas[, symbol := toupper(trimws(symbol))]
atlas <- atlas[!is.na(symbol) & symbol!=""][!duplicated(symbol)]
cat("atlas rows after dedup:", nrow(atlas), "\n")

rel <- fread(file.path(REL,"evidence_class_table.tsv"),
  select=c("symbol","primary_evidence_class","joint_testable"), showProgress=FALSE)
rel[, symbol := toupper(trimws(symbol))]
v <- merge(atlas, rel, by="symbol", all.x=TRUE)

frozen <- fread(file.path(REL,"evidence_class_validation_summary.tsv"), showProgress=FALSE)
lv <- c("neither","genetic_only","disease_state_only","convergent")

for (ep in c("proteomics","spatial_svg","single_cell")) {
  if (ep=="proteomics")  { tested <- !is.na(v$best_protein_padj); pos <- v$best_protein_padj < 0.05 }
  if (ep=="spatial_svg") { tested <- !is.na(v$spatial_is_svg);    pos <- v$spatial_is_svg %in% TRUE }
  if (ep=="single_cell") { tested <- !is.na(v$sc_best_padj);      pos <- v$sc_best_padj < 0.05 }
  keep <- tested & (v$joint_testable %in% TRUE)
  d <- v[keep]; d[, epos := pos[keep]]
  d <- d[primary_evidence_class %in% lv]
  s <- d[, .(recomputed_n_tested=.N, recomputed_n_positive=sum(epos,na.rm=TRUE)),
         by=primary_evidence_class]
  f <- frozen[endpoint==ep, .(primary_evidence_class, frozen_n_tested=n_tested,
                              frozen_n_positive=n_positive)]
  m <- merge(f, s, by="primary_evidence_class", all=TRUE)
  m <- m[match(lv, primary_evidence_class)]
  m[, match_tested   := frozen_n_tested   == recomputed_n_tested]
  m[, match_positive := frozen_n_positive == recomputed_n_positive]
  cat("\n===== ", ep, " =====\n"); print(m)
}
'
echo "=== DONE $(date) ==="
