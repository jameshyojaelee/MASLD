#!/usr/bin/env Rscript
# ─────────────────────────────────────────────────────────────────────────────
# Tier 1E (strengthening) — gene-set & cross-sectional "rewind" tests
#
# 1E single-gene signal is underpowered (n=15 regressors). Gene-set tests are
# more powerful. Two questions:
#   (A) Is regression a transcriptional REWIND of the established cross-sectional
#       progression axis?  -> correlate/GSEA regression-specific vs the canonical
#       disease-vs-control axis AND the per-stage fibrosis signatures.
#       Rewind  => strong NEGATIVE association. Distinct => ~0 / positive.
#   (B) What pathways define the regression-specific program? -> Hallmark fgsea.
#
# Inputs (RNA-seq/results/reversal/): de_regression_specific.csv,
#   de_progression_specific.csv, de_regress_change.csv, de_progress_change.csv
# Cross-sectional axes: integration/canonical_deg_results.csv (disease vs control)
#   + dream_results_stage_{steatosis,sh,cirrhosis}.csv (per-stage, if present)
# Outputs: rewind_test_summary.txt, fgsea_hallmark_{regression,progression}_specific.csv,
#   fgsea_progression_axis_in_regression.csv
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({ library(fgsea); library(msigdbr); library(data.table) })

root <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
rev  <- file.path(root, "RNA-seq/results/reversal")
intg <- file.path(root, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
rd   <- function(f) read.csv(file.path(rev, f), stringsAsFactors = FALSE)

rs <- rd("de_regression_specific.csv")
ps <- rd("de_progression_specific.csv")
canon <- read.csv(file.path(intg, "canonical_deg_results.csv"), stringsAsFactors = FALSE)

# ---- (A) rewind test vs cross-sectional disease axis (gene-level corr) -------
m <- merge(rs[, c("gene","symbol","logFC","t")],
           canon[, c("gene","logFC")], by = "gene", suffixes = c(".regSpec",".disease"))
mp <- merge(ps[, c("gene","logFC")], canon[, c("gene","logFC")], by = "gene",
            suffixes = c(".progSpec",".disease"))
rho_reg_dis  <- cor(m$logFC.regSpec,  m$logFC.disease,  method = "spearman", use="complete.obs")
rho_prog_dis <- cor(mp$logFC.progSpec, mp$logFC.disease, method = "spearman", use="complete.obs")

# ---- per-stage cross-sectional axes (if present) ----------------------------
stage_rho <- list()
for (s in c("steatosis","sh","cirrhosis")) {
  fp <- file.path(intg, sprintf("dream_results_stage_%s.csv", s))
  if (file.exists(fp)) {
    st <- read.csv(fp, stringsAsFactors = FALSE)
    lfc <- if ("logFC" %in% names(st)) "logFC" else grep("logFC|log2", names(st), value=TRUE)[1]
    a <- rs[, c("gene","logFC")]; names(a)[2] <- "reg"
    b <- st[, c("gene", lfc)];    names(b)[2] <- "stg"
    mm <- merge(a, b, by = "gene")
    stage_rho[[s]] <- cor(mm$reg, mm$stg, method="spearman", use="complete.obs")
  }
}

# ---- (A') GSEA: is the disease-UP / disease-DOWN set rewound in regression? --
canon_sig <- canon[!is.na(canon$padj) & canon$padj < 0.05 & abs(canon$logFC) > 0.5, ]
sets_axis <- list(
  disease_UP   = unique(canon_sig$symbol[canon_sig$logFC >  0]),
  disease_DOWN = unique(canon_sig$symbol[canon_sig$logFC <  0]))
rank_rs <- with(rs[!is.na(rs$t) & rs$symbol!="" & !is.na(rs$symbol), ],
                { v <- t; names(v) <- symbol; v[!duplicated(names(v))] })
rank_rs <- sort(rank_rs)
set.seed(42)
fg_axis <- fgsea(pathways = sets_axis, stats = rank_rs, eps = 0)
fwrite(fg_axis[, .(pathway, NES, pval, padj, size)],
       file.path(rev, "fgsea_progression_axis_in_regression.csv"))

# ---- (B) Hallmark fgsea on regression- and progression-specific rankings -----
hm <- msigdbr(species = "Homo sapiens", collection = "H")
hm_sets <- split(hm$gene_symbol, hm$gs_name)
rank_ps <- with(ps[!is.na(ps$t) & ps$symbol!="" & !is.na(ps$symbol), ],
                { v <- t; names(v) <- symbol; v[!duplicated(names(v))] })
rank_ps <- sort(rank_ps)
set.seed(42); fg_rs <- fgsea(hm_sets, rank_rs, eps = 0)
set.seed(42); fg_ps <- fgsea(hm_sets, rank_ps, eps = 0)
fwrite(fg_rs[order(padj), .(pathway,NES,pval,padj,size)], file.path(rev,"fgsea_hallmark_regression_specific.csv"))
fwrite(fg_ps[order(padj), .(pathway,NES,pval,padj,size)], file.path(rev,"fgsea_hallmark_progression_specific.csv"))

# ---- summary ----------------------------------------------------------------
sink(file.path(rev, "rewind_test_summary.txt"))
cat("Tier 1E strengthening — gene-set & rewind tests\n")
cat("================================================\n")
cat("(A) REWIND TEST vs cross-sectional axes (Spearman of logFC):\n")
cat(sprintf("    regression_specific  vs disease(canonical) : %+.3f\n", rho_reg_dis))
cat(sprintf("    progression_specific vs disease(canonical) : %+.3f  (sanity: should be +)\n", rho_prog_dis))
for (s in names(stage_rho)) cat(sprintf("    regression_specific  vs stage[%s]       : %+.3f\n", s, stage_rho[[s]]))
cat("    Rewind => regression_specific strongly NEGATIVE vs disease axis.\n")
cat("    Distinct resolution program => ~0 or positive.\n\n")
cat("(A') GSEA of disease-axis gene sets within the regression-specific ranking:\n")
print(fg_axis[, .(pathway, NES, padj, size)])
cat("    (rewind => disease_UP NES<0 & disease_DOWN NES>0)\n\n")
cat("(B) Top Hallmark pathways — regression_specific (padj<0.25):\n")
print(fg_rs[padj < 0.25][order(padj), .(pathway, NES, padj, size)])
cat("\n(B) Top Hallmark pathways — progression_specific (padj<0.25):\n")
print(fg_ps[padj < 0.25][order(padj), .(pathway, NES, padj, size)])
sink()
cat(readLines(file.path(rev, "rewind_test_summary.txt")), sep = "\n")
