#!/usr/bin/env Rscript
# benchmark_cutoffs.R
# ---------------------------------------------------------------------------
# Benchmark filtering cutoffs for the Cas13 library, on the canonical STAR -s 2
# human dream mega-analysis. CONTROL-RECOVERY-CENTRIC, UP-DIRECTIONAL (library
# selects UP-in-disease genes for a knockdown screen).
#
# Significance is FIXED at 0.05; the EFFECT-SIZE cutoff is the variable (x-axis,
# discrete points). Two conditions:
#   RAW  : padj < 0.05 & log2FC      > x
#   ASHR : lfsr < 0.05 & shrunk_log2FC > x
# (x=0 = no magnitude floor = pure padj<0.05 / pure lfsr<0.05.)
# NOTE raw log2FC and ashr shrunk are different scales; same x is not equivalent.
#
# Metrics (per condition x cutoff): recall + fold-enrichment of the ~1,025-gene
# OpenTargets MASLD gold standard (independent), mouse cross-species precision,
# and 5-fold leave-one-cohort-out (LOCO) reproducibility. n_genes for size context.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
has_pw <- requireNamespace("patchwork", quietly = TRUE)
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INTDIR <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
FIGDIR <- file.path(BASE, "Cas13_Library_Design/figures")
DATDIR <- file.path(BASE, "Cas13_Library_Design/data")
strip_v <- function(x) sub("[.][0-9]+$", "", x)

# ---- 1. dream + per-gene flags ---------------------------------------------
d <- fread(file.path(INTDIR, "dream_results_ashr.csv"))
if (!"se" %in% names(d)) d[, se := abs(logFC / t)]
d[, gb := strip_v(gene)]
d[, full_sign := sign(logFC)]

# OpenTargets MASLD gold standard (~1,025 genes; aggregated genetics/drugs/literature
# -> INDEPENDENT of our transcriptomics). Used for recall + fold-enrichment (panels A/B).
ot <- fread(file.path(BASE, "data/published_gene_panels/opentargets_masld_2025.tsv"),
            skip = "gene_symbol")
ot_genes <- unique(ot$gene_symbol)
d[, is_ot := symbol %in% ot_genes]
N_OT    <- sum(d$is_ot)            # OpenTargets genes present in the dream universe
OT_BASE <- N_OT / nrow(d)          # genome-wide prevalence (enrichment denominator)
cat(sprintf("OpenTargets MASLD genes in universe: %d / %d (baseline %.1f%%)\n", N_OT, nrow(d), 100 * OT_BASE))

atlas <- fread(file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"),
               select = c("ensembl_id", "mouse_meta_padj"))
atlas[, gb := strip_v(ensembl_id)]
d <- merge(d, atlas[, .(gb, mouse_padj = mouse_meta_padj)], by = "gb", all.x = TRUE)
d[, is_mouse := !is.na(mouse_padj) & mouse_padj < 0.05]

loo_files <- list.files(file.path(INTDIR, "loo_cv"),
                        pattern = "^dream_loo_GSE.*\\.csv$", full.names = TRUE)
lm <- NULL
for (f in loo_files) {
  fo <- fread(f); gcol <- intersect(c("gene","Gene","row"), names(fo))[1]
  lcol <- intersect(c("logFC","log2FoldChange"), names(fo))[1]
  pcol <- intersect(c("padj","adj.P.Val"), names(fo))[1]
  fo[, gb := strip_v(get(gcol))]
  m <- merge(d[, .(gb, full_sign)], fo[, .(gb, lf = get(lcol), pa = get(pcol))], by = "gb", all.x = TRUE)
  v <- as.integer(!is.na(m$lf) & sign(m$lf) == m$full_sign & !is.na(m$pa) & m$pa < 0.05)
  lm <- if (is.null(lm)) v else lm + v
}
d[, loco_repro := lm / length(loo_files)]

# ---- 2. evaluate the two conditions across a grid of effect-size cutoffs ----
metrics_of <- function(sel) { s <- d[sel]; list(
  n_genes    = nrow(s),
  recall_ot  = sum(s$is_ot) / N_OT,        # sensitivity: fraction of OpenTargets MASLD genes recovered
  prec_ot    = mean(s$is_ot),              # precision: fraction of selected that are OpenTargets
  enrich_ot  = mean(s$is_ot) / OT_BASE,    # fold-enrichment vs genome baseline
  prec_mouse = mean(s$is_mouse),
  loco       = mean(s$loco_repro, na.rm = TRUE)) }

CUTS <- c(0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.75, 1.0)  # 0.15 added (candidate cutoff)
RAW  <- "padj<0.05 & log2FC > x"
ASHR <- "lfsr<0.05 & shrunk > x"
res <- rbindlist(c(
  lapply(CUTS, function(x) c(list(condition = RAW,  cutoff = x),
            metrics_of(d[, !is.na(padj) & padj < 0.05 & logFC > x]))),
  lapply(CUTS, function(x) c(list(condition = ASHR, cutoff = x),
            metrics_of(d[, !is.na(lfsr) & lfsr < 0.05 & shrunk_logFC > x])))))
res[, condition := factor(condition, levels = c(RAW, ASHR))]
pal <- setNames(c("#B2182B", "#3B4CC0"), c(RAW, ASHR))
fwrite(res, file.path(DATDIR, "cutoff_benchmark_by_effectsize.csv"))
cat("\n=== metrics by effect-size cutoff (significance fixed at 0.05) ===\n"); print(res)

# ---- 3. figures: metric vs effect-size cutoff, two conditions --------------
bt <- theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(),
        legend.title = element_blank(), legend.position = "top",
        plot.title = element_text(face = "bold", size = 11))
pct <- scales::percent_format(accuracy = 1)
xlab <- "Effect-size cutoff"
mk <- function(y, ytitle, title) ggplot(res, aes(cutoff, get(y), color = condition)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2.2) +
  scale_color_manual(values = pal) + scale_x_continuous(breaks = CUTS) +
  scale_y_continuous(labels = pct) + labs(x = xlab, y = ytitle, title = title) + bt

pA <- mk("recall_ot", "Recall of OpenTargets MASLD genes",
         sprintf("A  Recall of OpenTargets MASLD genes (gold standard, n=%d)", N_OT))

pB <- ggplot(res, aes(cutoff, enrich_ot, color = condition)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50", linewidth = 0.4) +
  geom_line(linewidth = 0.8) + geom_point(size = 2.2) +
  scale_color_manual(values = pal) + scale_x_continuous(breaks = CUTS) +
  labs(x = xlab, y = "Fold-enrichment vs genome",
       title = sprintf("B  Enrichment for OpenTargets MASLD genes (baseline %.1f%%)", 100 * OT_BASE)) + bt

pC <- mk("prec_mouse", "Fraction mouse-replicated", "C  Mouse cross-species precision")
pD <- mk("loco", "Mean LOCO reproducibility", "D  Cross-cohort reproducibility (5-fold LOCO)")

panels <- list(A = pA, B = pB, C = pC, D = pD)
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)
for (nm in names(panels))
  ggsave(file.path(FIGDIR, sprintf("S_lib_8_cutoff_benchmark_%s.pdf", nm)),
         panels[[nm]], width = 6.8, height = 4.2, useDingbats = FALSE)
cat("\nWrote individual panels A-D to", FIGDIR, "\n")
