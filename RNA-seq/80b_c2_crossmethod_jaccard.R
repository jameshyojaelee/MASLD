#!/usr/bin/env Rscript
# 80b_c2_crossmethod_jaccard.R
#
# C2 RECOUNT add-on (2026-06-08): cross-method concordance of the MuSiC
# hepatocyte-intrinsic attribution between the OLD dream-canonical run and the
# NEW C2 (limma-voom-qw) recount.
#
# The historical "BayesPrism-MuSiC Jaccard 0.774" compared two *deconvolution
# algorithms* (BayesPrism vs MuSiC) on the SAME bulk DEG set. A global
# disease-vs-control BayesPrism DE no longer exists on disk (only F-stage
# transition BayesPrism runs remain), so that exact comparison cannot be
# recomputed cheaply. What IS cheap and directly relevant to this task is the
# DEG-method-swap concordance: does the hepatocyte-intrinsic gene set move when
# the bulk DEG spine changes dream -> C2? That is what we compute here.
#
# Env: rnaseq
# Outputs: RNA-seq/results/celltype_attribution/c2_recount/

suppressPackageStartupMessages({ library(data.table) })

BASE   <- Sys.getenv("MASLD_PROJECT_ROOT",
                     "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
CAUSAL <- file.path(BASE, "RNA-seq/results/causal_inference")
OUTDIR <- file.path(BASE, "RNA-seq/results/celltype_attribution/c2_recount")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

jaccard <- function(a, b) {
  a <- unique(a); b <- unique(b)
  length(intersect(a, b)) / length(union(a, b))
}
base_ensg <- function(x) sub("\\.\\d+$", "", x)

# ---- OLD dream MuSiC attribution (current canonical on disk) ----
old <- fread(file.path(CAUSAL, "deconv_attribution_scores.csv"))
old[, ensg := base_ensg(gene)]
old_hep <- old[category == "Hepatocyte_intrinsic", unique(ensg)]

# ---- NEW C2 MuSiC attribution (this recount) ----
new <- fread(file.path(CAUSAL, "c2_recount", "deconv_attribution_scores.csv"))
new[, ensg := base_ensg(gene)]
new_hep <- new[category == "Hepatocyte_intrinsic", unique(ensg)]

# ---- Concordance ----
inter <- length(intersect(old_hep, new_hep))
uni   <- length(union(old_hep, new_hep))
jac   <- jaccard(old_hep, new_hep)

# Category-level cross-tab on the shared gene universe
m <- merge(old[, .(ensg, cat_dream = category)],
           new[, .(ensg, cat_c2    = category)],
           by = "ensg", all = FALSE)
xtab <- dcast(m, cat_dream ~ cat_c2, value.var = "ensg", fun.aggregate = length)

lines <- c(
  "=== C2 deconvolution recount: hepatocyte-intrinsic cross-method concordance ===",
  "",
  sprintf("OLD dream Hepatocyte_intrinsic (deconv_attribution_scores.csv): %d genes", length(old_hep)),
  sprintf("NEW C2    Hepatocyte_intrinsic (c2_recount/...):                %d genes", length(new_hep)),
  sprintf("Intersection: %d", inter),
  sprintf("Union:        %d", uni),
  sprintf("Jaccard(dream-hep-intrinsic, C2-hep-intrinsic) = %.3f", jac),
  "",
  "Category cross-tab (rows = dream category, cols = C2 category):",
  capture.output(print(xtab)),
  "",
  "NOTE: The historical 'BayesPrism-MuSiC Jaccard 0.774' compared two",
  "      deconvolution algorithms on one DEG set; a global disease-vs-control",
  "      BayesPrism DE no longer exists on disk (only F-stage BayesPrism runs),",
  "      so that exact number is not recomputable here. The Jaccard above is the",
  "      DEG-method-swap concordance (dream -> C2) of the MuSiC hep-intrinsic set."
)
writeLines(lines, file.path(OUTDIR, "crossmethod_jaccard_summary.txt"))
writeLines(lines)

fwrite(data.table(
  metric = c("old_dream_hep_intrinsic", "new_c2_hep_intrinsic",
             "intersection", "union", "jaccard"),
  value  = c(length(old_hep), length(new_hep), inter, uni, round(jac, 4))
), file.path(OUTDIR, "crossmethod_jaccard.csv"))

cat("\nDone. Wrote crossmethod_jaccard_summary.txt + crossmethod_jaccard.csv\n")
