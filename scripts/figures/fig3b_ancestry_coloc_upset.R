#!/usr/bin/env Rscript
# Standalone renderer for Fig 3 Panel B (cross-ancestry COLOC UpSet).
# UpSet plot of eGene intersections across {EUR, EAS, AFR, SAS} at PP4>0.5,
# best SuSiE PP.H4 per gene per ancestry (ABF fallback when SuSiE did not
# converge). Renders without UpSetR/ComplexUpset (pure ggplot+patchwork).

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- file.path(FIG3_DIR, "panels")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)

EUR_17 <- c(
  "2019_31311600_NAFLD_EUR", "2020_32298765_NAFLD_EUR",
  "2021_34128465_PDFF_EUR",  "2021_34841290_NAFLD_EUR",
  "2021_34957434_PDFF_EUR",  "2022_36402844_PDFF_EUR",
  "2023_36280732_NAFLD_deCode_EUR", "2023_36280732_NAFLD_Intermountain_EUR",
  "2023_36280732_NAFLD_UKBB_EUR", "FinnGen_NAFLD", "FinnGen_NASH",
  "UKBB_ALT", "UKBB_AST", "UKBB_GGT"
)
BBJ_5 <- c("BBJ_ALT", "BBJ_AST", "BBJ_GGT")
AFR_3 <- c("PanUKBB_AFR_ALT", "PanUKBB_AFR_AST", "PanUKBB_AFR_GGT")
SAS_3 <- c("PanUKBB_CSA_ALT", "PanUKBB_CSA_AST", "PanUKBB_CSA_GGT")

ancestry_for_gwas <- function(g) {
  fcase(
    g %in% EUR_17, "EUR",
    g %in% BBJ_5,  "EAS",
    g %in% AFR_3,  "AFR",
    g %in% SAS_3,  "SAS",
    default       = NA_character_
  )
}

PP4_THRESH <- 0.5

# 1KG-everywhere assembly (LD-panel-fair):
#   EUR (14 GWAS): susie_coloc_1kg/susie_coloc_all_gwas_1kg.csv (EUR re-run on
#     1KG after UKBB-sghatan retired 2026-04-23). UKBB↔1KG concordance r²=0.998.
#   EAS / AFR / SAS: susie_coloc/susie_coloc_all_gwas.csv (already 1KG EAS/AFR/SAS).
sc_eur <- fread(file.path(BASE,
  "GWAS/finemapping/results/susie_coloc_1kg/susie_coloc_all_gwas_1kg.csv"))
sc_other <- fread(file.path(BASE,
  "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv"))

sc_eur   <- sc_eur[gwas_name %in% EUR_17]
sc_other <- sc_other[gwas_name %in% c(BBJ_5, AFR_3, SAS_3)]

# Align columns (sc_eur and sc_other share PP.H4.abf / PP.H4.susie); keep both.
common_cols <- intersect(names(sc_eur), names(sc_other))
sc <- rbindlist(list(sc_eur[, ..common_cols], sc_other[, ..common_cols]),
                use.names = TRUE)

sc[, ancestry := ancestry_for_gwas(gwas_name)]
sc <- sc[!is.na(ancestry)]
sc[, pp4_best := pmax(PP.H4.susie, PP.H4.abf, na.rm = TRUE)]
sc[is.infinite(pp4_best), pp4_best := NA_real_]
cat(sprintf("[fig3b] LD panel: 1KG everywhere (EUR=%d rows from 1kg re-run, non-EUR=%d rows from primary file)\n",
            nrow(sc_eur), nrow(sc_other)))

# Best PP4 per gene per ancestry across that ancestry's GWAS portfolio
per_gene_anc <- sc[!is.na(pp4_best),
                   .(max_pp4 = max(pp4_best, na.rm = TRUE)),
                   by = .(gene, ancestry)]

ancestry_levels <- c("EUR", "EAS", "AFR", "SAS")
sig <- per_gene_anc[max_pp4 > PP4_THRESH]

# Build gene × ancestry membership table
gene_sets <- dcast(sig, gene ~ ancestry, value.var = "max_pp4",
                   fill = NA_real_)
for (a in ancestry_levels) {
  if (!a %in% names(gene_sets)) gene_sets[, (a) := NA_real_]
  gene_sets[, (paste0("in_", a)) := !is.na(get(a))]
}

# Encode each gene's intersection (binary string)
membership_cols <- paste0("in_", ancestry_levels)
gene_sets[, intersection_id := apply(.SD, 1,
  function(r) paste0(as.integer(r), collapse = "")),
  .SDcols = membership_cols]

intersections <- gene_sets[, .(n_genes = .N), by = intersection_id]
# Decode bitmask back to set membership
intersections[, c(membership_cols) := lapply(seq_along(membership_cols), function(i) {
  substr(intersection_id, i, i) == "1"
})]
intersections[, n_sets := rowSums(.SD), .SDcols = membership_cols]
intersections <- intersections[n_sets > 0]
setorder(intersections, -n_genes)
intersections[, ix_rank := .I]
intersections[, ix_label := factor(ix_rank, levels = ix_rank)]

# Per-ancestry totals (set sizes)
set_sizes <- data.table(
  ancestry = ancestry_levels,
  n_genes  = sapply(ancestry_levels, function(a) sum(gene_sets[[paste0("in_", a)]]))
)
set_sizes[, ancestry := factor(ancestry, levels = rev(ancestry_levels))]

# --- Top panel: intersection bar -------------------------------------------
ancestry_colors <- c(
  "EUR" = "#C9265E", "EAS" = "#F4511E",
  "AFR" = "#00695C", "SAS" = "#7B1FA2"
)
intersections[, primary := fcase(
  n_sets >= 4, "All 4",
  n_sets == 3, "3-way",
  n_sets == 2, "2-way",
  n_sets == 1 & in_EUR, "EUR only",
  n_sets == 1 & in_EAS, "EAS only",
  n_sets == 1 & in_AFR, "AFR only",
  n_sets == 1 & in_SAS, "SAS only"
)]
primary_colors <- c(
  "All 4"     = "#212121",
  "3-way"     = "#455A64",
  "2-way"     = "#90A4AE",
  "EUR only"  = ancestry_colors["EUR"],
  "EAS only"  = ancestry_colors["EAS"],
  "AFR only"  = ancestry_colors["AFR"],
  "SAS only"  = ancestry_colors["SAS"]
)

p_top <- ggplot(intersections, aes(x = ix_label, y = n_genes, fill = primary)) +
  geom_col(width = 0.75, color = NA) +
  geom_text(aes(label = n_genes), vjust = -0.3, size = GEOM_TEXT_6PT, color = "gray20") +
  scale_fill_manual(values = primary_colors, name = NULL,
                    breaks = c("All 4", "3-way", "2-way",
                               "EUR only", "EAS only", "AFR only", "SAS only")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = NULL, y = "eGenes (intersection size)") +
  theme_masld() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right",
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6),
        plot.margin = margin(2, 4, 1, 4))

# --- Bottom-right: dot matrix ----------------------------------------------
dot_dt <- melt(intersections[, c("ix_label", membership_cols), with = FALSE],
               id.vars = "ix_label",
               variable.name = "ancestry", value.name = "in_set")
dot_dt[, ancestry := factor(sub("^in_", "", ancestry),
                            levels = rev(ancestry_levels))]

# Connector segments: for each intersection with >=2 sets, span min..max
segs <- dot_dt[in_set == TRUE,
               .(ymin = min(as.integer(ancestry)),
                 ymax = max(as.integer(ancestry))),
               by = ix_label]
segs <- segs[ymin != ymax]

p_dots <- ggplot(dot_dt, aes(x = ix_label, y = ancestry)) +
  geom_point(aes(color = in_set), size = 2.0) +
  { if (nrow(segs)) geom_segment(data = segs,
      aes(x = ix_label, xend = ix_label,
          y = ymin, yend = ymax),
      inherit.aes = FALSE,
      color = "#212121", linewidth = 0.5) } +
  scale_color_manual(values = c("TRUE" = "#212121", "FALSE" = "#E0E0E0"),
                     guide = "none") +
  labs(x = NULL, y = NULL) +
  theme_masld() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(face = "plain", size = 6),
        plot.margin = margin(1, 4, 4, 4))

# --- Bottom-left: per-ancestry set-size bar (rotated) ----------------------
p_left <- ggplot(set_sizes, aes(y = ancestry, x = n_genes, fill = as.character(ancestry))) +
  geom_col(width = 0.7, color = NA) +
  geom_text(aes(label = n_genes), hjust = 1.1, size = GEOM_TEXT_6PT, color = "white") +
  scale_fill_manual(values = ancestry_colors, guide = "none") +
  scale_x_reverse(expand = expansion(mult = c(0.05, 0))) +
  labs(x = "Set size", y = NULL) +
  theme_masld() +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(1, 1, 4, 4))

# --- Layout: top spans full width; bottom = (left | dots) ------------------
spacer <- plot_spacer()

bottom <- p_left + p_dots + plot_layout(widths = c(1, 5))
top    <- spacer + p_top + plot_layout(widths = c(1, 5))

message(sprintf(
  "[caption] Cross-ancestry COLOC eGene intersections (PP.H4 > %.1f). 23-GWAS portfolio: 14 EUR + 3 EAS + 3 AFR + 3 SAS Pan-UKBB; best SuSiE PP.H4 per gene (ABF fallback).",
  PP4_THRESH))

upset <- top / bottom +
  plot_layout(heights = c(2.4, 1.0))

out_pdf <- file.path(PANEL_DIR, "ancestry_coloc_upset.pdf")
out_csv <- file.path(FIG3_DIR, "ancestry_coloc_upset.csv")
save_fig(upset, out_pdf, width = fig_full_width * 0.55, height = 3.2)

# Persist intersection table for reviewers / caption
fwrite(intersections[, .(intersection_id, in_EUR, in_EAS, in_AFR, in_SAS,
                          n_sets, n_genes, primary)],
       out_csv)
cat("[fig3b] Wrote:", out_pdf, "\n")
cat("[fig3b] Wrote:", out_csv, "\n")
print(intersections[, .(in_EUR, in_EAS, in_AFR, in_SAS, n_genes, primary)])
cat("\nSet sizes:\n")
print(set_sizes)
