#!/usr/bin/env Rscript
# ============================================================================
# deg_meta_heterogeneity_ridge.R
# SUPP (integration_value) — cross-cohort heterogeneity of integrated DEGs
#
# Density of metafor meta_I2 for Tier-1 integrated DEGs vs the rest of the
# tested transcriptome. Integrated DEGs are markedly LESS heterogeneous across
# cohorts (median I^2 ~50%) than the transcriptome at large (~73%). High-I2
# anchors (FASN/SCD/FGF21) are only estimable BECAUSE of integration.
#
# FRAME: meta_I2 comes from the RETAINED metafor random-effects SENSITIVITY
# arm, NOT the C2 canonical limma-voom-qw DEG-calling model. It is used here
# purely as a cross-cohort consistency diagnostic for genes called by C2.
#
# Output: <FIGS_INTVAL_DIR>/deg_meta_heterogeneity_ridge.pdf  (84 x 64 mm)
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

dir.create(FIGS_INTVAL_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_PDF <- file.path(FIGS_INTVAL_DIR, "deg_meta_heterogeneity_ridge.pdf")
OUT_CSV <- file.path(FIGS_INTVAL_DIR, "deg_meta_heterogeneity_ridge.csv")

INT <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")

# ---------------------------------------------------------------------------
# Load metafor sensitivity-arm heterogeneity + C2 canonical DEG calls
# ---------------------------------------------------------------------------
meta <- fread(file.path(INT, "meta_analysis_results.csv"))
deg  <- fread(file.path(INT, "canonical_deg_results.csv"))

# I^2 only defined where >=2 cohorts contribute
meta <- meta[!is.na(meta_I2) & n_datasets >= 2]
# meta gene IDs are UNVERSIONED ENSG; canonical/anchors are VERSIONED -> strip
meta[, gene_base := sub("\\..*$", "", gene)]

# Tier-1 integrated DEG (TREAT canonical): treat_fdr<0.05 at lfc=0.25
deg[, gene_base := sub("\\..*$", "", gene)]
deg[, is_deg := is_dream_deg(deg)]
deg_genes <- deg[is_deg == TRUE, gene_base]

meta[, deg_status := ifelse(gene_base %in% deg_genes,
                            "Integrated DEG", "Other tested genes")]
meta[, deg_status := factor(deg_status,
                            levels = c("Other tested genes", "Integrated DEG"))]

# ---------------------------------------------------------------------------
# Medians + hero numbers
# ---------------------------------------------------------------------------
med <- meta[, .(med_I2 = median(meta_I2)), by = deg_status]
med_deg   <- med[deg_status == "Integrated DEG",     med_I2]
med_other <- med[deg_status == "Other tested genes", med_I2]
cat(sprintf("[hero] median meta_I2: Integrated DEGs=%.1f%% (n=%d) vs Other=%.1f%% (n=%d)\n",
            med_deg,   meta[deg_status == "Integrated DEG", .N],
            med_other, meta[deg_status == "Other tested genes", .N]))

# ---------------------------------------------------------------------------
# High-I2 anchor genes (need integration to estimate at all)
# ---------------------------------------------------------------------------
anchors <- c("FASN", "SCD", "FGF21")
anchor_ens <- sub("\\..*$", "", deg[symbol %in% anchors, gene])
anch <- meta[gene_base %in% anchor_ens]
# attach symbol for labelling
anch <- merge(anch, unique(deg[symbol %in% anchors,
              .(gene_base = sub("\\..*$", "", gene), symbol)]),
              by = "gene_base", all.x = TRUE)
anch[, gene := symbol]
anch[, y_lab := 0]   # placed on the rug baseline
cat("[hero] anchors: ",
    paste(sprintf("%s I2=%.1f%%", anch$gene, anch$meta_I2), collapse = "; "), "\n")

fwrite(meta[, .(gene, meta_I2, n_datasets, meta_logFC, deg_status)], OUT_CSV)

# ---------------------------------------------------------------------------
# Colors: Liang — DEGs magenta, background neutral gray
# ---------------------------------------------------------------------------
fill_cols <- c("Integrated DEG" = "#C9265E", "Other tested genes" = "#9E9E9E")

p <- ggplot(meta, aes(x = meta_I2, fill = deg_status, color = deg_status)) +
  geom_density(alpha = 0.35, linewidth = 0.5, adjust = 1.1) +
  geom_vline(xintercept = med_deg,   color = "#C9265E",
             linetype = "dashed", linewidth = 0.45) +
  geom_vline(xintercept = med_other, color = "#616161",
             linetype = "dashed", linewidth = 0.45) +
  # high-I2 anchors as rug + labels
  geom_rug(data = anch, aes(x = meta_I2), inherit.aes = FALSE,
           color = "#1565C0", linewidth = 0.6, length = unit(0.06, "npc")) +
  geom_text_repel(data = anch, aes(x = meta_I2, y = 0.004, label = gene),
                  inherit.aes = FALSE, size = GEOM_TEXT_6PT, color = "#1565C0",
                  fontface = "italic", segment.size = 0.25,
                  direction = "y", nudge_y = 0.006, min.segment.length = 0) +
  annotate("text", x = med_deg, y = Inf, label = sprintf("DEG median\n%.1f%%", med_deg),
           hjust = 1.1, vjust = 1.4, size = GEOM_TEXT_6PT, color = "#C9265E", lineheight = 0.9) +
  annotate("text", x = med_other, y = Inf, label = sprintf("genome-wide\n%.1f%%", med_other),
           hjust = -0.08, vjust = 1.4, size = GEOM_TEXT_6PT, color = "#616161", lineheight = 0.9) +
  scale_fill_manual(values = fill_cols, name = NULL) +
  scale_color_manual(values = fill_cols, name = NULL) +
  scale_x_continuous(name = expression("Cross-cohort heterogeneity  " * I^2 * "  (%)"),
                     limits = c(0, 100), expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(name = "Density", expand = expansion(mult = c(0, 0.06))) +
  theme_masld(base_size = 6) +
  theme(
    legend.position = c(0.99, 0.62),
    legend.justification = c(1, 0.5),
    legend.text     = element_text(size = 6),
    legend.key.size = unit(0.16, "cm"),
    legend.background = element_blank()
  )

ggsave(OUT_PDF, p, width = 84 / 25.4, height = 64 / 25.4,
       units = "in", device = cairo_pdf)
cat(sprintf("[saved] %s\n", OUT_PDF))
message("[caption] Integrated DEGs are less heterogeneous across cohorts (metafor sensitivity-arm I^2, not the C2 DEG model) — a consistency diagnostic")
