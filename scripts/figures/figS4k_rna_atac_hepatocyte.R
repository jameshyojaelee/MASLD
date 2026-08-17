#!/usr/bin/env Rscript
# ==============================================================================
# figS4k — hepatocyte mRNA↔promoter-accessibility TREND (snATAC companion to 4b).
# Same scatter grammar as Fig 4b (mRNA↔protein), but y = donor-level hepatocyte
# promoter-accessibility log2FC (snATAC, GSE244832). Each point is a gene; x =
# mRNA log2FC (pooled bulk RNA-seq, C2).
#
# THIS IS A DIRECTIONAL TREND, NOT A SIGNIFICANCE CALL — deliberately. The
# n=18-donor multiome (5 NORMAL controls) is underpowered for per-gene DA: 0
# genes pass BH FDR<0.10, and the nominal-p hit-rate among prioritized promoters
# (5.5%) is indistinguishable from the 5% null (1.1×). So we make NO per-gene
# significance claim on the accessibility axis. What survives at n=18 is the
# ACTIVITY TREND (the paper's own donor-level verdict: "trends survive,
# significance collapses"): mRNA and promoter accessibility are directionally
# coupled, and that coupling is ENRICHED in the prioritized universe over the
# measured background (Spearman rho 0.35 vs 0.24; 64% vs 58% direction-concordant,
# both above the 50% chance level). That enrichment — not any single gene — is the
# claim.
#
# RECONCILIATION WITH FIG 4d: 4d deliberately EXCLUDES this disease-DA contrast
# from its GENETICS analysis (variant→caQTL direction: the 5-control disease-DA is
# caQTL-anti-concordant). Here the disease-DA is used ONLY as an mRNA-direction
# trend — no genetic-direction and no per-gene claim — so that caQTL liability
# does not transfer; the "underpowered, 5 controls" caveat still does, which is
# exactly why no significance gate is applied.
#
#   scatter: x = mRNA log2FC (bulk RNA-seq); y = hepatocyte promoter-accessibility
#            log2FC (snATAC donor pseudobulk, GSE244832; gene-level = the
#            TSS-proximal called Hep peak within 2 kb of the TSS). light-gray =
#            measured background; teal = prioritized universe (fig4a/4d). Dashed
#            fits over each group make the enriched (steeper) prioritized slope
#            visible; the prioritized-vs-background rho gap is annotated in-panel.
#
# Source: Analysis/ATAC/Human_Multiome/results/l8_annotated_corrected/
#           scatac_da_gene_annotated.csv  (donor-pseudobulk Hep DA, 07b)
#         canonical_deg_results.csv  (bulk RNA-seq logFC, C2)
#         prioritized_universe_FINAL.txt  (fig4a/4d shared target set)
# Output: figures/main/fig5_molecular_context/panels/figS4k.pdf   (PDF only)
# Env:    rnaseq
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
set.seed(42)
FAM  <- "Helvetica"
BASE <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PROMOTER_WINDOW <- 2000L   # bp; standard promoter definition for gene-level accessibility

# ── ATAC: gene-level hepatocyte promoter-accessibility log2FC ──────────────────
# The annotated DA keeps peaks out to their nearest-gene distance (median ~10 kb);
# restrict to TSS-proximal promoter peaks and collapse to the closest one per gene.
atac <- fread(file.path(BASE, "Analysis/ATAC/Human_Multiome/results/l8_annotated_corrected/scatac_da_gene_annotated.csv"))
atac <- atac[distance_to_tss <= PROMOTER_WINDOW & is.finite(logFC)]
setorder(atac, gene_symbol, distance_to_tss)
prom <- atac[, .SD[1L], by = gene_symbol][, .(gene = gene_symbol, atac_logFC = logFC)]

# ── RNA: bulk mRNA log2FC (canonical C2) ──────────────────────────────────────
rna <- fread(file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv"))
rna <- unique(rna[is.finite(logFC), .(gene = symbol, rna_logFC = logFC)], by = "gene")

con <- merge(prom, rna, by = "gene")
con <- con[is.finite(atac_logFC) & is.finite(rna_logFC)]

uni <- trimws(readLines(file.path(BASE, "Analysis/Spatial/results/universe_validation/prioritized_universe_FINAL.txt")))
uni <- uni[uni != ""]
con[, in_uni := gene %in% uni]

# ── trend statistics: prioritized vs background (Spearman + direction concordance) ─
sgn_conc <- function(a, b) mean(sign(a) == sign(b))
du <- con[in_uni == TRUE]         # prioritized measured set
db <- con[in_uni == FALSE]        # measured background
rho_uni  <- cor(du$rna_logFC, du$atac_logFC, method = "spearman"); conc_uni <- sgn_conc(du$rna_logFC, du$atac_logFC)
rho_bg   <- cor(db$rna_logFC, db$atac_logFC, method = "spearman"); conc_bg  <- sgn_conc(db$rna_logFC, db$atac_logFC)
n_uni <- nrow(du); n_bg <- nrow(db)
cat(sprintf("[figS4k] measured=%d | background n=%d rho=%.3f conc=%.1f%% | prioritized n=%d rho=%.3f conc=%.1f%%\n",
            nrow(con), n_bg, rho_bg, 100*conc_bg, n_uni, rho_uni, 100*conc_uni))

# colours / point styling (two groups only — NO significance/direction partition)
COL <- c(background = "#E0E0E0", prioritized = "#00695C")
con[, grp := factor(fifelse(in_uni, "prioritized", "background"),
                    levels = c("background", "prioritized"))]
setorder(con, grp)   # draw background first, prioritized on top

# ── scatter (single panel; enrichment carried by the two fit lines + annotation) ─
ann <- sprintf("ρ = %.2f prioritized · %.2f background\n%.0f%% direction-concordant",
               rho_uni, rho_bg, 100*conc_uni)
xr <- range(con$rna_logFC); yr <- range(con$atac_logFC)
sc <- ggplot(con, aes(rna_logFC, atac_logFC, colour = grp, alpha = grp)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.22, colour = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.22, colour = "grey70") +
  # both clouds kept light (they overlap heavily near 0,0); the trend lives in the fits
  rasterize_layer(geom_point(size = 0.4, stroke = 0), dpi = 450) +
  # OLS fits over each group: the steeper prioritized slope IS the enrichment
  geom_smooth(data = db, aes(rna_logFC, atac_logFC), method = "lm", formula = y ~ x, se = FALSE,
              linetype = "dashed", linewidth = 0.5, colour = "grey45", inherit.aes = FALSE) +
  geom_smooth(data = du, aes(rna_logFC, atac_logFC), method = "lm", formula = y ~ x, se = FALSE,
              linewidth = 0.7, colour = "#00695C", inherit.aes = FALSE) +
  scale_colour_manual(values = COL, labels = c(background = "Measured background", prioritized = "Prioritized target"),
                      name = NULL, guide = guide_legend(override.aes = list(size = 1.4, alpha = 1), ncol = 1)) +
  scale_alpha_manual(values = c(background = 0.30, prioritized = 0.38), guide = "none") +
  annotate("text", x = xr[1], y = yr[2], hjust = 0, vjust = 1, size = 6/.pt, family = FAM,
           lineheight = 0.9, label = ann) +
  labs(x = expression("mRNA "*log[2]*"FC (bulk RNA-seq)"),
       y = expression("hepatocyte promoter "*log[2]*"FC (snATAC, GSE244832)")) +
  theme_masld() +
  theme(legend.position = c(0.99, 0.02), legend.justification = c(1, 0),
        legend.background = element_rect(fill = alpha("white", 0.7), colour = NA),
        legend.key.size = unit(0.28, "lines"), legend.text = element_text(size = 6),
        legend.margin = margin(1, 2, 1, 2), plot.margin = margin(2, 2, 1, 2))

# ── output (single scatter; no marginal) ───────────────────────────────────────
out <- file.path(FIG4_DIR, "panels", "figS4k.pdf")
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
ggsave(out, sc, width = 3.0, height = 2.13, device = grDevices::cairo_pdf)
cat("[figS4k] saved:", out, "\n")

message(sprintf(paste0(
  "CAPTION (Fig S4k): Hepatocyte mRNA↔promoter-accessibility TREND — the snATAC companion to Fig 4b. Each point ",
  "is a gene; x = mRNA log2FC (pooled bulk RNA-seq, C2), y = hepatocyte promoter-accessibility log2FC (snATAC ",
  "donor pseudobulk, GSE244832; gene-level = the TSS-proximal called hepatocyte peak within %d bp of the TSS, ",
  "per-donor OLS on log2(CPM+1), n=18 donors incl. 5 NORMAL controls). Light-gray = measured background ",
  "(n=%s); teal = the prioritized universe shared across Fig 4 (n=%s measured of the transcriptomic 8,088 union ",
  "genetic 3,038 target set). mRNA and promoter accessibility are directionally coupled, and the coupling is ",
  "ENRICHED in the prioritized set over background: Spearman rho=%.2f vs %.2f, and %.0f%% vs %.0f%% of genes ",
  "agree in direction (both above the 50%% chance level; dashed lines are per-group OLS fits, the steeper ",
  "prioritized slope shown at right). IMPORTANT — this is a DIRECTIONAL TREND, NOT a significance call: the ",
  "18-donor multiome is underpowered for per-gene differential accessibility (0 genes pass BH FDR<0.10; the ",
  "nominal-p hit-rate among prioritized promoters, 5.5%%, is ~1.1x the 5%% null), so NO per-gene accessibility ",
  "claim is made — consistent with the paper's donor-level verdict that activity TRENDS survive while per-gene ",
  "significance collapses. Fig 4d excludes this disease-DA contrast from its GENETICS analysis (the 5-control ",
  "disease-DA is caQTL-anti-concordant); here it is used ONLY as an mRNA-direction trend, so that liability does ",
  "not transfer, while the underpowering caveat is exactly why no significance gate is applied. Static ",
  "accessibility breadth and the genetic-colocalization view are in Fig 4d; sequence-level motif disruption in ",
  "Fig S4j."),
  PROMOTER_WINDOW, format(n_bg, big.mark=","), format(n_uni, big.mark=","),
  rho_uni, rho_bg, 100*conc_uni, 100*conc_bg))
