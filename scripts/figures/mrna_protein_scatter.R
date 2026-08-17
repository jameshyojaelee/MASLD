#!/usr/bin/env Rscript
# KEY MESSAGE (fig4A): atlas DEGs replicate at the protein level in BOTH liver
# tissue and plasma, two compartments that do not always agree. Per gene, mRNA
# log2FC (canonical limma-voom quality-weighted C2) vs protein log2FC (DIA-MS,
# MASLD vs control). Two panels: liver tissue (PXD051911, Boel et al.) and plasma
# (PXD052937, Sourianarayanane et al.). Concordance is stronger in tissue than in
# plasma, and a handful of genes (e.g. COL5A1) flip sign between compartments.
#
# Both panels are built from one source, with DEG / conserved flags derived the
# same way, and the plasma "all-genes" rho reproduces the canonical 0.31 (text);
# liver reproduces 0.34. Source: protein_transcript_concordance_v3.csv.
#
# Output: figures/main/fig5_molecular_context/mrna_protein_scatter.pdf
# Env:    rnaseq

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(ggrepel)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# ── Gene-level flags (dataset-independent), derived once for both panels ──────
# DEG set = canonical bulk Tier-1 at the proteomics-concordance threshold
# (padj < 0.05, |logFC| > 0.3); Conserved = cross-species conserved core (>=3).
deg_tab <- read.csv(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv"),
  stringsAsFactors = FALSE)
deg_set <- unique(deg_tab$symbol[!is.na(deg_tab$padj) & deg_tab$padj < 0.05 &
                                 abs(deg_tab$logFC) > 0.3])
cons <- read.csv(file.path(BASE,
  "Analysis/Cross_Species_Concordance/results/concordance_atlas_unified.csv"),
  stringsAsFactors = FALSE)
cons_set <- unique(cons$human_symbol[suppressWarnings(as.numeric(cons$n_concordant)) >= 3])

conc <- read.csv(file.path(BASE,
  "Analysis/Proteomics/results/protein_transcript_concordance_v3.csv"),
  stringsAsFactors = FALSE)

CASE <- c("CYP3A4", "HKDC1")
rho <- function(x, y) suppressWarnings(cor(x, y, method = "spearman"))

build_panel <- function(ds_id, ptitle) {
  d <- conc %>%
    filter(dataset == ds_id) %>%
    transmute(gene, mrna_lfc = bulk_logFC, prot_lfc = protein_logFC,
              prot_padj = protein_padj) %>%
    filter(!is.na(mrna_lfc), !is.na(prot_lfc)) %>%
    mutate(is_deg   = gene %in% deg_set,
           is_cons  = gene %in% cons_set,
           prot_sig = !is.na(prot_padj) & prot_padj < 0.05,
           is_both  = is_deg & prot_sig,
           dirmatch = sign(mrna_lfc) == sign(prot_lfc))

  r_all  <- rho(d$mrna_lfc, d$prot_lfc)
  r_deg  <- rho(d$mrna_lfc[d$is_deg],  d$prot_lfc[d$is_deg])
  r_cons <- rho(d$mrna_lfc[d$is_cons], d$prot_lfc[d$is_cons])
  nb <- sum(d$is_both); nm <- sum(d$is_both & d$dirmatch)
  pct <- if (nb) 100 * nm / nb else NA
  message(sprintf("[%s] all=%.2f (n=%d) | DEG=%.2f (n=%d) | cons=%.2f (n=%d) | both %d/%d (%.0f%%)",
                  ds_id, r_all, nrow(d), r_deg, sum(d$is_deg), r_cons, sum(d$is_cons), nm, nb, pct))

  bg     <- filter(d, !is_deg)          # non-DEG background (gray)
  deg    <- filter(d,  is_deg)          # DEGs (magenta), including the both-significant ones
  case_d <- filter(d, gene %in% CASE)   # the two named case genes (ringed + labeled)

  bg_layer <- geom_point(data = bg, aes(mrna_lfc, prot_lfc),
                         color = masld_colors[["ns"]], size = 0.4, alpha = 0.28)
  if (requireNamespace("ggrastr", quietly = TRUE))
    bg_layer <- ggrastr::rasterise(bg_layer, dpi = 300)

  xr <- range(d$mrna_lfc); yr <- range(d$prot_lfc)

  ggplot() +
    geom_hline(yintercept = 0, linewidth = 0.25, color = "gray80") +
    geom_vline(xintercept = 0, linewidth = 0.25, color = "gray80") +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.3, linetype = "dashed", color = "gray60") +
    bg_layer +
    geom_point(data = deg, aes(mrna_lfc, prot_lfc),
               color = masld_colors[["up"]], size = 0.8, alpha = 0.6) +
    # Label ONLY the two named case genes, with a leader line/arrow to the point
    # (no ring). No both-sig colour, legend, or stat box; rho values are in the
    # caption. Gene-label size (geom units) ~= theme size 13pt / 2.845.
    geom_text_repel(data = case_d, aes(mrna_lfc, prot_lfc, label = gene),
                    fontface = "bold.italic", size = 13 / 2.845,
                    color = "gray15", segment.color = "gray30", segment.size = 0.4,
                    min.segment.length = 0, box.padding = 1.1, point.padding = 0.4,
                    arrow = arrow(length = unit(0.02, "npc"), type = "closed"),
                    max.overlaps = Inf, seed = 1) +
    labs(x = "mRNA log2FC", y = "Protein log2FC", title = ptitle) +
    coord_cartesian(xlim = c(xr[1] - 0.05, xr[2] + 0.05),
                    ylim = c(yr[1] - 0.1, yr[2] + 0.1)) +
    theme_masld() + theme_pub()
}

pL <- build_panel("PXD051911", "Liver tissue")
pR <- build_panel("PXD052937", "Plasma")

# All visible text at ONE uniform size (13pt): overall title, panel titles,
# axis titles, tick labels (and the gene labels, set in geom units above).
SZ <- 13
p_out <- (pL | pR) +
  plot_annotation(
    title = "Atlas DEGs replicate at the protein level in liver and plasma",
    theme = theme(plot.title = element_text(size = SZ, face = "plain"))
  ) &
  theme(legend.position = "none",
        plot.title  = element_text(size = SZ, face = "plain"),
        axis.title  = element_text(size = SZ),
        axis.text   = element_text(size = SZ, color = "black"))

out <- file.path(FIG4_DIR, "_supp", "mrna_protein_scatter.pdf")
cairo_pdf(out, width = fig_full_width, height = fig_full_width * 0.6)
print(p_out)
invisible(dev.off())
message("Saved: ", out)
