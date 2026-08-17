#!/usr/bin/env Rscript
# fig3_deg_upset_cohorts.R
# UpSet of DEG sets: our INTEGRATED pooled human signature (limma-voom quality-
# weighted C2) vs the 5 control-bearing cohorts' OWN per-study DE (the same
# voomWithQualityWeights limma-voom, 02_per_study_de.R) — plus the union /
# replication structure of the integrated signature across cohorts.
#
# All 6 sets thresholded identically on the CANONICAL TREAT scale:
# TREAT FDR < 0.05 at lfc = 0.25 — the new canonical Tier-1 definition
# (= 1,918 pooled DEGs; 1,419 up / 499 down). treat() tests H0:|log2FC|<=0.25,
# so the effect floor is folded INTO the test (no separate |logFC| cut). The
# integrated set carries treat_fdr from canonical_deg_results.csv; each per-study
# set is given an analytical TREAT here (limma::treat reconstructed from logFC +
# SE + df.total, identical to the engine in rebuild_cas13_library.R) and gated
# identically.
# Output: figures/main/fig3_RNAseq/panels/  (main fig3 panels; moved here 2026-06-11)
#   figs3_deg_upset.pdf            UpSet of integrated + 5 per-study sets
#   figs3b_deg_replication.pdf      integrated DEGs by # cohorts replicating
#   deg_upset_summary.csv          set sizes + union overlaps

suppressPackageStartupMessages({
  library(data.table)
  library(UpSetR)
  library(ggplot2)
  library(grid)
  library(ashr)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results")
PER  <- file.path(INT, "per_study")
OUTDIR <- file.path(BASE, "figures/main/fig3_RNAseq/panels")   # main fig3 panels dir
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
suppressWarnings(tryCatch(source(file.path(BASE, "scripts/figures/publication_theme.R")),
                          error = function(e) NULL))

# Canonical DEG definition (2026-08-12): padj < 0.05 AND |log2FC| > 0.5.
# BOTH arms of this UpSet use the same padj cut and the same effect-size floor.
# Do not let them drift apart: if the per-cohort arms are held to a stricter floor
# than the integrated arm, every gene between the two floors is barred from all
# per-cohort sets while being admitted to the integrated set, and can therefore
# only land in "Integrated only" — manufacturing the panel's headline.
CANON_PADJ_CUT <- 0.05
CANON_LFC_CUT  <- 0.50
strip_ver <- function(x) sub("\\.[0-9]+$", "", as.character(x))

# Per-cohort canonical DEGs: returns gene + effect SIGN under the same gate the
# integrated arm uses. (Superseded the per-cohort analytical TREAT reconstruction
# on 2026-08-12, when TREAT stopped being the canonical definition.)
canonical_degs <- function(d) {
  p <- intersect(c("padj", "adj.P.Val"), names(d))[1]
  if (is.na(p)) stop("canonical_degs: no adjusted p-value column")
  sig <- d[!is.na(get(p)) & !is.na(logFC) &
             get(p) < CANON_PADJ_CUT & abs(logFC) > CANON_LFC_CUT]
  unique(sig[, .(gene = strip_ver(gene), sign = sign(logFC))], by = "gene")
}

# 5 control-bearing cohorts pooled into the integrated disease-vs-control signature
# (UpSet set labels displayed by accession)
cohorts <- data.table(
  gse   = c("GSE126848","GSE130970","GSE135251","GSE162694","GSE213621"),
  label = c("GSE126848","GSE130970","GSE135251","GSE162694","GSE213621"))

# ── Integrated pooled limma-voom-qw C2 (canonical treat_fdr) ─────────────────
can <- fread(file.path(INT, "integration/canonical_deg_results.csv"))
can_dt <- unique(can[!is.na(padj) & !is.na(logFC) &
                       padj < CANON_PADJ_CUT & abs(logFC) > CANON_LFC_CUT,
                     .(gene = strip_ver(gene), sign = sign(logFC))], by = "gene")
can_deg <- can_dt$gene
can_sign <- setNames(can_dt$sign, can_dt$gene)
cat(sprintf("Integrated (limma-voom-qw C2): %d DEG at padj<%.2f & |log2FC|>%.2f\n",
            length(can_deg), CANON_PADJ_CUT, CANON_LFC_CUT))

# ── Per-study limma-voom-qw, SAME canonical gate as the integrated arm ───────
sets <- list(Integrated = can_deg)
cohort_sign <- list()   # gene -> effect sign, per cohort (for direction-matched support)
for (i in seq_len(nrow(cohorts))) {
  f <- file.path(PER, paste0(cohorts$gse[i], "_de_results.csv"))
  d <- fread(f)
  sig <- canonical_degs(d)
  sets[[cohorts$label[i]]] <- sig$gene
  cohort_sign[[cohorts$label[i]]] <- setNames(sig$sign, sig$gene)
  cat(sprintf("  %-8s (%s): %d DEG (analytical TREAT)\n",
              cohorts$label[i], cohorts$gse[i], length(sig$gene)))
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
final_pdf <- file.path(OUTDIR, "figs3_deg_upset.pdf")
raw_pdf   <- file.path(OUTDIR, ".deg_upset_raw.pdf")
pdf(raw_pdf, width = fig_full_width, height = fig_full_width * 6.5 / 11, useDingbats = FALSE)
upset(mat,
      sets = rev(c("Integrated", cohorts$label)),
      keep.order = TRUE,
      order.by = "freq", nintersects = 30,
      mainbar.y.label = "DEGs in intersection",
      sets.x.label = "DEGs per set (FDR<0.05, lfc=0.25)",
      main.bar.color = "#D41159", sets.bar.color = "#1A85FF",
      matrix.color = "#333333", shade.color = "#E8E8E8",
      text.scale = c(1.6, 1.4, 1.4, 1.2, 1.5, 1.3),
      point.size = 2.8, line.size = 0.9,
      mb.ratio = c(0.62, 0.38))
dev.off()
# UpSetR's internal grid.newpage() leaves a blank leading page; drop it (keep page 2+).
gs <- Sys.which("gs")
if (nzchar(gs) && system2(gs, c("-sDEVICE=pdfwrite","-dFirstPage=2","-dNOPAUSE","-dBATCH","-dQUIET",
        paste0("-sOutputFile=", shQuote(final_pdf)), shQuote(raw_pdf))) == 0 && file.exists(final_pdf)) {
  file.remove(raw_pdf)
} else file.rename(raw_pdf, final_pdf)
cat(sprintf("\nWrote %s\n", final_pdf))
message("[caption] DEG concordance: integrated signature vs each cohort's own limma-voom-qw DE")

# ── Panel: integrated DEGs by # of cohorts replicating (compact bar) ─────────
# Direction-matched cohort support on the CANONICAL TREAT scale: for each
# integrated Tier-1 DEG (TREAT FDR<0.05, lfc=0.25), a cohort "supports" it if
# that cohort's within-cohort analytical TREAT (same engine) calls it at the
# same gate AND in the same direction as the integrated effect.
support_n <- sapply(C, function(g) {
  sg <- can_sign[g]
  sum(vapply(cohorts$label, function(l) {
    s <- cohort_sign[[l]][g]   # NA if gene absent in this cohort's DEG set
    isTRUE(!is.na(s) && s == sg)
  }, logical(1)))
})
rd <- data.table(n_cohorts = 0:5)[
  data.table(n_cohorts = factor(support_n, levels = 0:5))[, .N, by = n_cohorts][
    , .(n_cohorts = as.integer(as.character(n_cohorts)), n_genes = N)],
  on = "n_cohorts"]
rd[is.na(n_genes), n_genes := 0]
setorder(rd, n_cohorts)
rd[, pct := 100 * n_genes / sum(n_genes)]
stopifnot(sum(rd$n_genes) == length(C))   # all integrated Tier-1 DEGs accounted for
pal <- colorRampPalette(c("#9E9E9E", "#1A85FF", "#D41159"))(6)
pRep <- ggplot(rd, aes(factor(n_cohorts), n_genes, fill = factor(n_cohorts))) +
  geom_col(width = 0.82) +
  geom_text(aes(label = sprintf("%d\n(%.0f%%)", n_genes, pct)),
            vjust = -0.2, size = GEOM_TEXT_6PT, lineheight = 0.82, fontface = "plain") +
  scale_fill_manual(values = pal, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +
  labs(x = "# of 5 cohorts replicating", y = "Integrated DEGs") +
  theme_bw(base_size = 13) +
  theme(axis.title         = element_text(size = 6),
        axis.text          = element_text(size = 6),
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        plot.margin        = margin(4, 6, 4, 4))
ggsave(file.path(OUTDIR, "figs3b_deg_replication.pdf"), pRep,
       width = 4.6, height = 3.2, useDingbats = FALSE)
cat(sprintf("Wrote %s\n", file.path(OUTDIR, "figs3b_deg_replication.pdf")))
message(sprintf("[caption] Integrated DEGs (n=%d) replicating per cohort", length(C)))
cat("\n=== done ===\n")
