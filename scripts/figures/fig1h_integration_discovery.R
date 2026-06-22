#!/usr/bin/env Rscript
# fig1h_integration_discovery.R — Per-gene integrated logFC vs. cohort concordance count.
# Highlights genes that reach significance only in the dream mega-analysis ("integration-only"),
# vs. genes already significant in 2+ individual cohorts ("shared discovery").
#
# x: number of cohorts (0..5) where padj<0.05 AND sign(logFC) == sign(integrated logFC)
# y: integrated dream logFC
# Cohort universe matches Fig 1e UpSet: Suppli, Hoang, Govaere, Bril, Chen.
# Output: figures/main/fig3_RNAseq/panels/fig3a_integration_discovery.pdf

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(ggforce)   # geom_sina
  library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- file.path(FIG2_DIR, "panels")   # relocated fig1 -> fig3_RNAseq (mirrors fig1h_optionb.R)
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)

PADJ_INT       <- 0.05   # integrated significance threshold (panel inclusion)
# Deliberately NO |logFC| cutoff on the integrated side: integration's power is
# detecting small-effect genes single cohorts miss; an LFC cutoff preferentially
# filters out the integration-only population.
PADJ_COHORT    <- 0.05   # per-cohort padj for x-axis "concordantly significant"
LFC_CANONICAL  <- 0.5    # per-cohort |logFC| for canonical DEG status (Integration-only color)
N_LABEL_GENES   <- 4    # top integration-only genes by padj
N_BIN0_OUTLIERS <- 2    # additional bin-0 genes by largest |logFC|

FIVE_COHORTS <- c("GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621")

# ----------------------------------------------------------------------------
# Load data and compute per-gene cohort-concordance count
# ----------------------------------------------------------------------------
message("Loading data...")
dream <- load_dream_results()
ps    <- load_per_study_de()[dataset %in% FIVE_COHORTS]

# Strip Ensembl version for robust joining (per-study and dream may differ)
dream[, gene_clean := sub("\\..*", "", gene)]
ps[,    gene_clean := sub("\\..*", "", gene)]

dream_lookup <- dream[, .(gene_clean, bulk_logFC, bulk_padj)]

ps_join <- merge(ps[, .(gene_clean, dataset, logFC, padj)],
                 dream_lookup, by = "gene_clean", all.x = FALSE)
# x-axis: cohorts with concordant nominal significance (padj + same direction)
ps_join[, concordant_sig := !is.na(padj) & padj < PADJ_COHORT &
                            !is.na(logFC) & !is.na(bulk_logFC) &
                            sign(logFC) == sign(bulk_logFC) & sign(logFC) != 0]
# color: cohorts that meet canonical DEG threshold (padj + |logFC|>0.5 + direction)
ps_join[, canonical_deg := concordant_sig & abs(logFC) > LFC_CANONICAL]

concordance <- ps_join[, .(
  n_cohorts_concordant = sum(concordant_sig, na.rm = TRUE),
  n_cohorts_canonical  = sum(canonical_deg, na.rm = TRUE)
), by = gene_clean]

panel_df <- merge(dream, concordance, by = "gene_clean", all.x = TRUE)
panel_df[is.na(n_cohorts_concordant), n_cohorts_concordant := 0L]
panel_df[is.na(n_cohorts_canonical),  n_cohorts_canonical := 0L]

# Restrict the panel to integrated DEGs (statistical significance only).
# Non-DEGs are not part of the question. See the LFC_INT note above for why we
# do NOT impose |logFC|>0.5 here.
panel_df <- panel_df[!is.na(bulk_logFC) & !is.na(bulk_padj) &
                     bulk_padj < PADJ_INT]

# Integration-only = no cohort meets canonical DEG criteria for this gene
# (padj < 0.05 AND |logFC| > 0.5 AND same direction). Genes can still be
# Integration-only at high x-axis values if cohorts agree on direction at
# padj < 0.05 but none had |logFC| > 0.5 — orthogonal to the x-axis.
panel_df[, category := fifelse(n_cohorts_canonical == 0,
                               "Integration-only", "Canonical DEG in ≥1 cohort")]

# ----------------------------------------------------------------------------
# Summary stats
# ----------------------------------------------------------------------------
cat_counts <- panel_df[, .N, by = category]
setorder(cat_counts, -N)
message("Category counts:")
for (i in seq_len(nrow(cat_counts))) {
  message(sprintf("  %-20s %s", cat_counts$category[i], comma(cat_counts$N[i])))
}

bin_counts <- panel_df[, .N, by = n_cohorts_concordant]
setorder(bin_counts, n_cohorts_concordant)
message("\nPer-bin DEG counts (n_cohorts_concordant: n DEGs):")
for (i in seq_len(nrow(bin_counts))) {
  message(sprintf("  bin %d: %s DEGs",
                  bin_counts$n_cohorts_concordant[i],
                  comma(bin_counts$N[i])))
}

n_int_only   <- cat_counts[category == "Integration-only", N]
n_replicated <- cat_counts[category == "Canonical DEG in ≥1 cohort", N]
if (length(n_int_only)   == 0) n_int_only   <- 0
if (length(n_replicated) == 0) n_replicated <- 0

# ----------------------------------------------------------------------------
# Per-gene supplementary table
# ----------------------------------------------------------------------------
out_table <- panel_df[, .(gene, symbol, bulk_logFC, bulk_padj,
                          n_cohorts_concordant, n_cohorts_canonical, category)]
setorder(out_table, bulk_padj)
fwrite(out_table, file.path(PANEL_DIR, "fig1h_concordance_table.csv"))
message(sprintf("\nWrote table: %s rows", comma(nrow(out_table))))

# ----------------------------------------------------------------------------
# Plot
# ----------------------------------------------------------------------------
message("\nGenerating plot...")

cat_colors <- c(
  "Integration-only"     = "#C9265E",   # Liang deep magenta — the spotlight
  "Canonical DEG in ≥1 cohort" = "#9E9E9E"    # neutral gray — contextual background, not a competing category
)

panel_df[, x_factor := factor(n_cohorts_concordant, levels = 0:5)]
panel_df[, abs_logFC := abs(bulk_logFC)]
panel_df[, direction := fifelse(bulk_logFC >= 0, "Up", "Down")]
panel_df[, category := factor(category,
                              levels = c("Canonical DEG in ≥1 cohort", "Integration-only"))]
setorder(panel_df, category)   # so Integration-only renders on top

# Curated labels: disease-relevant integration-only genes spanning the spectrum
# of MASLD biology (drug targets, canonical hepatocyte regulators, fibrosis
# ECM, antioxidant/detox, immune). GSTM1 and HLA-DQB1 are also the two bin-0
# |logFC| outliers, so they double as the "high-effect integration-only" callouts.
curated_labels <- c(
  "GSTM1",      # bin 0, |logFC|=1.4, antioxidant/detox; biggest integration-only effect
  "HLA-DQB1",   # bin 0, |logFC|=1.1, immune/MHC class II
  "NR1H4",      # bin 1, FXR — obeticholic acid target
  "THRB"        # bin 2, thyroid hormone receptor β — resmetirom target
)
label_df <- panel_df[symbol %in% curated_labels &
                     category == "Integration-only"]
setorder(label_df, bulk_padj)

# Per-bin n label
n_label_df <- panel_df[, .(N = .N),
                       by = .(x_factor, n_cohorts_concordant)]
setorder(n_label_df, n_cohorts_concordant)
y_top <- max(panel_df$abs_logFC, na.rm = TRUE) * 1.06

p <- ggplot(panel_df,
            aes(x = x_factor, y = abs_logFC,
                color = category, group = x_factor)) +
  # Sina: shared density per bin — `group = x_factor` overrides the implicit
  # grouping by `color`, so Integration-only and Per-study-replicated dots
  # share one merged distribution shape with magenta layered on top.
  rasterize_layer(
    geom_sina(maxwidth = 0.85, size = 0.5, alpha = 0.55, shape = 16,
              seed = 42)
  ) +
  # Slim boxplot for summary stats per bin
  geom_boxplot(width = 0.12, fill = scales::alpha("white", 0.85),
               color = "gray30", outlier.shape = NA, linewidth = 0.3,
               show.legend = FALSE) +
  # Top integration-only gene labels (use abs_logFC for y-position)
  geom_text_repel(data = label_df,
                  aes(x = x_factor, y = abs(bulk_logFC), label = symbol),
                  inherit.aes = FALSE,
                  size = 2.8, color = "black",
                  fontface = "italic",
                  segment.size = 0.2, segment.color = "gray50",
                  box.padding = 0.35, point.padding = 0.15,
                  force = 3, max.overlaps = Inf,
                  min.segment.length = 0, seed = 42,
                  show.legend = FALSE) +
  # Reference line at |logFC| = 0.5 — canonical Tier 1 DEG cutoff.
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.3, color = "gray60") +
  scale_color_manual(values = cat_colors,
                     breaks = c("Integration-only", "Canonical DEG in ≥1 cohort"),
                     name = NULL) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.16)),
                     limits = c(0, NA)) +
  labs(
    x = "Cohorts with concordant significance (padj < 0.05, of 5)",
    y = expression("Integrated |log"[2]*" FC|"),
    title = "Effect size scales with cross-cohort replication"
  ) +
  theme_masld(base_size = 10) +
  theme(
    plot.title    = element_text(size = 11, face = "bold", margin = margin(b = 1)),
    plot.subtitle = element_text(size = 8.5, color = "gray30", margin = margin(b = 2)),
    plot.margin   = margin(2, 3, 2, 2),
    axis.title.x  = element_text(margin = margin(t = 1)),
    axis.title.y  = element_text(margin = margin(r = 1)),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.background = element_blank(),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.spacing = unit(2, "pt"),
    legend.key.height = unit(0.3, "cm")
  )

save_fig(p, file.path(PANEL_DIR, "fig3a_integration_discovery.pdf"),
         width = fig_half_width * 1.55, height = 3.5)

fp <- file.path(PANEL_DIR, "fig3a_integration_discovery.pdf")
if (file.exists(fp)) {
  message(sprintf("\nOutput: %s (%s)", fp,
                  utils:::format.object_size(file.size(fp), "auto")))
} else {
  message("ERROR: Output file not created")
}
