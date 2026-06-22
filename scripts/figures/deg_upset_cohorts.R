#!/usr/bin/env Rscript
# fig3_deg_upset_cohorts.R
# UpSet of DEG sets: our INTEGRATED pooled human signature (limma-voom quality-
# weighted C2) vs the 5 control-bearing cohorts' OWN per-study DE (the same
# voomWithQualityWeights limma-voom, 02_per_study_de.R) — plus the union /
# replication structure of the integrated signature across cohorts.
#
# All 6 sets thresholded identically: padj < 0.05 AND |logFC| > 0.5 — the same
# Tier-1 raw cutoff used to define the canonical/integrated DEGs (NUMBERS.md).
# (padj alone, no LFC floor, gives 13,043 integrated genes; the |logFC|>0.5
#  floor brings it to the canonical 1,853.)
# Output: figures/main/fig3_RNAseq/panels/  (main fig3 panels; moved here 2026-06-11)
#   fig3b_deg_upset.pdf            UpSet of integrated + 5 per-study sets
#   figs3b_deg_replication.pdf      integrated DEGs by # cohorts replicating
#   deg_upset_summary.csv          set sizes + union overlaps

suppressPackageStartupMessages({
  library(data.table)
  library(UpSetR)
  library(ggplot2)
  library(grid)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results")
PER  <- file.path(INT, "per_study")
OUTDIR <- file.path(BASE, "figures/main/fig3_RNAseq/panels")   # main fig3 panels dir
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
suppressWarnings(tryCatch(source(file.path(BASE, "scripts/figures/publication_theme.R")),
                          error = function(e) NULL))

PADJ <- 0.05
LFC  <- 0.5   # match canonical Tier-1 raw DEG definition: padj<0.05 & |logFC|>0.5
strip_ver <- function(x) sub("\\.[0-9]+$", "", as.character(x))

# 5 control-bearing cohorts pooled into the integrated disease-vs-control signature
# (UpSet set labels displayed by accession)
cohorts <- data.table(
  gse   = c("GSE126848","GSE130970","GSE135251","GSE162694","GSE213621"),
  label = c("GSE126848","GSE130970","GSE135251","GSE162694","GSE213621"))

# ── Integrated pooled limma-voom-qw C2 ──────────────────────────────────────
can <- fread(file.path(INT, "integration/canonical_deg_results.csv"))
can_deg <- unique(strip_ver(can[!is.na(padj) & !is.na(logFC) &
                                padj < PADJ & abs(logFC) > LFC, gene]))
cat(sprintf("Integrated (limma-voom-qw C2): %d DEG at padj<%.2f & |logFC|>%.1f\n",
            length(can_deg), PADJ, LFC))

# ── Per-study limma-voom-qw (same method + cutoff, 02_per_study_de.R) ────────
sets <- list(Integrated = can_deg)
for (i in seq_len(nrow(cohorts))) {
  f <- file.path(PER, paste0(cohorts$gse[i], "_de_results.csv"))
  d <- fread(f)
  deg <- unique(strip_ver(d[!is.na(adj.P.Val) & !is.na(logFC) &
                            adj.P.Val < PADJ & abs(logFC) > LFC, gene]))
  sets[[cohorts$label[i]]] <- deg
  cat(sprintf("  %-8s (%s): %d DEG\n", cohorts$label[i], cohorts$gse[i], length(deg)))
}

# ── Union / overlap summary ─────────────────────────────────────────────────
per_union <- unique(unlist(sets[cohorts$label]))
C <- sets$Integrated
summ <- data.table(
  metric = c("integrated_DEG","per_study_union(any cohort)","integrated AND union",
             "integrated only (no cohort replicates)","union only (not in integrated)",
             paste0("integrated replicated in >=", 1:5, " cohorts")),
  n = c(length(C), length(per_union), length(intersect(C, per_union)),
        length(setdiff(C, per_union)), length(setdiff(per_union, C)),
        sapply(1:5, function(k) sum(rowSums(sapply(cohorts$label, function(l) C %in% sets[[l]])) >= k))))
fwrite(summ, file.path(OUTDIR, "deg_upset_summary.csv"))
cat("\nUnion/overlap summary:\n"); print(summ)

# ── Panel: UpSet (integrated + 5 per-study) ─────────────────────────────────
mat <- fromList(sets)
final_pdf <- file.path(OUTDIR, "fig3b_deg_upset.pdf")
raw_pdf   <- file.path(OUTDIR, ".deg_upset_raw.pdf")
pdf(raw_pdf, width = 11, height = 6.5, useDingbats = FALSE)
upset(mat,
      sets = rev(c("Integrated", cohorts$label)),
      keep.order = TRUE,
      order.by = "freq", nintersects = 30,
      mainbar.y.label = "DEGs in intersection",
      sets.x.label = "DEGs per set (padj<0.05, |logFC|>0.5)",
      main.bar.color = "#D41159", sets.bar.color = "#1A85FF",
      matrix.color = "#333333", shade.color = "#E8E8E8",
      text.scale = c(1.6, 1.4, 1.4, 1.2, 1.5, 1.3),
      point.size = 2.8, line.size = 0.9,
      mb.ratio = c(0.62, 0.38))
grid.text("DEG concordance: integrated signature vs each cohort's own limma-voom-qw DE",
          x = 0.66, y = 0.98, gp = gpar(fontsize = 12, fontface = "bold"))
dev.off()
# UpSetR's internal grid.newpage() leaves a blank leading page; drop it (keep page 2+).
gs <- Sys.which("gs")
if (nzchar(gs) && system2(gs, c("-sDEVICE=pdfwrite","-dFirstPage=2","-dNOPAUSE","-dBATCH","-dQUIET",
        paste0("-sOutputFile=", shQuote(final_pdf)), shQuote(raw_pdf))) == 0 && file.exists(final_pdf)) {
  file.remove(raw_pdf)
} else file.rename(raw_pdf, final_pdf)
cat(sprintf("\nWrote %s\n", final_pdf))

# ── Panel: integrated DEGs by # of cohorts replicating (compact bar) ─────────
# Cohort support uses the CANONICAL A2 definition the manuscript is locked to:
# for each integrated Tier-1 DEG, a cohort "supports" it if that cohort's
# within-cohort LVQW re-fit (same engine as the integrated analysis) calls it at
# padj<0.05 & |logFC|>0.5 AND in the same direction as the integrated effect.
# Source of truth: recompute_integration_only_cohort_support_C2.R ->
# integration_only_C2/A1_cohort_support_breakdown_C2.csv  (0 cohorts = 91).
# The previous version counted the standalone per_study DE CSVs direction-agnostic
# (0 cohorts = 44, i.e. 98% supported), which disagreed with the manuscript's
# 1,762 / 95.1% — a different per-cohort producer, not a different finding.
A1 <- fread(file.path(INT, "integration", "integration_only_C2",
                      "A1_cohort_support_breakdown_C2.csv"))
rd <- A1[, .(n_cohorts = n_cohorts_supporting, n_genes)]
rd[, pct := 100 * n_genes / sum(n_genes)]
stopifnot(sum(rd$n_genes) == length(C))   # 1,853 integrated Tier-1 DEGs
pal <- colorRampPalette(c("#9E9E9E", "#1A85FF", "#D41159"))(6)
pRep <- ggplot(rd, aes(factor(n_cohorts), n_genes, fill = factor(n_cohorts))) +
  geom_col(width = 0.82) +
  geom_text(aes(label = sprintf("%d\n(%.0f%%)", n_genes, pct)),
            vjust = -0.2, size = 3.4, lineheight = 0.82, fontface = "bold") +
  scale_fill_manual(values = pal, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +
  labs(title = sprintf("Integrated DEGs (n=%d) replicating per cohort", length(C)),
       x = "# of 5 cohorts replicating", y = "Integrated DEGs") +
  theme_bw(base_size = 13) +
  theme(plot.title         = element_text(face = "bold", size = 12),
        axis.title         = element_text(size = 12),
        axis.text          = element_text(size = 11),
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        plot.margin        = margin(4, 6, 4, 4))
ggsave(file.path(OUTDIR, "figs3b_deg_replication.pdf"), pRep,
       width = 4.6, height = 3.2, useDingbats = FALSE)
cat(sprintf("Wrote %s\n", file.path(OUTDIR, "figs3b_deg_replication.pdf")))
cat("\n=== done ===\n")
