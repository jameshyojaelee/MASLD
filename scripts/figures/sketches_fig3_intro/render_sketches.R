#!/usr/bin/env Rscript
# Renders 10 layout options + 4 hybrids for a Fig 3 intro/summary panel.
# Each design saved as a separate PDF in figures/sketches/fig3_intro/.
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggrepel); library(patchwork)
  library(ggalluvial); library(ggforce); library(scales)
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
OUT  <- file.path(BASE, "figures/sketches/fig3_intro")
DAT  <- file.path(BASE, "scripts/figures/sketches_fig3_intro/data")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
source(file.path(BASE, "scripts/figures/publication_theme.R"))

# ---- Load tables ----
portfolio  <- fread(file.path(DAT, "gwas_portfolio.csv"))
cascade    <- fread(file.path(DAT, "method_cascade.csv"))
method_tot <- fread(file.path(DAT, "method_totals.csv"))
eqtl       <- fread(file.path(DAT, "eqtl_panels.csv"))
ancestry   <- fread(file.path(DAT, "ancestry_summary.csv"))

portfolio[, ancestry := factor(ancestry, levels = c("EUR","EAS","AFR","SAS"))]
portfolio[, trait := factor(trait, levels = c(
  "NAFLD/NASH","HCC","Cirrhosis","PDFF (MRI)","Liver enzymes","Other"))]
ancestry[, ancestry := factor(ancestry, levels = c("EUR","EAS","AFR","SAS"))]

ancestry_colors <- c(EUR = "#1F77B4", EAS = "#D62728", AFR = "#2CA02C", SAS = "#9467BD")
trait_colors    <- c("NAFLD/NASH"   = "#8C564B",
                     "HCC"          = "#E377C2",
                     "Cirrhosis"    = "#BCBD22",
                     "PDFF (MRI)"   = "#17BECF",
                     "Liver enzymes"= "#7F7F7F",
                     "Other"        = "#AEC7E8")
method_colors   <- c("SuSiE (single ancestry)" = "#377EB8",
                     "SuSiE-X (cross-ancestry)" = "#E41A1C",
                     "meSuSiE (multi-ethnic)"   = "#4DAF4A")

save_sk <- function(p, name, w = 8, h = 5) {
  ggsave(file.path(OUT, name), p, width = w, height = h, device = cairo_pdf)
  cat(" wrote", name, "\n")
}

cat("=== Rendering Fig 3 intro sketch options ===\n")

# =========================================================================
# Option 1: Portfolio Dashboard (3-panel horizontal)
# =========================================================================
p1a <- ggplot(portfolio, aes(x = N_tot, y = reorder(label, N_tot))) +
  geom_segment(aes(xend = 0, yend = label, color = ancestry), linewidth = 0.4) +
  geom_point(aes(color = ancestry, size = pmax(susie_n_loci_converged, 1))) +
  scale_x_log10(labels = scales::label_comma()) +
  scale_size_continuous(range = c(1.5, 5), name = "SuSiE loci") +
  scale_color_manual(values = ancestry_colors) +
  labs(x = "Sample size (log)", y = NULL, title = "GWAS portfolio (28 studies)") +
  theme_masld() +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 6))

p1b <- ggplot(eqtl, aes(x = panel, y = n_samples, fill = panel)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = paste0("N=", scales::label_comma()(n_samples), "\n",
                               scales::label_comma()(n_egenes), " eGenes")),
            vjust = -0.4, size = 2.5) +
  scale_fill_manual(values = c("Broadaway liver" = "#FFA726",
                               "GTEx v8 Liver"   = "#FB8C00")) +
  labs(x = NULL, y = "Samples", title = "eQTL panels") +
  theme_masld() + theme(legend.position = "none") +
  expand_limits(y = 1500)

p1c <- ggplot(method_tot, aes(x = method, y = n_genes_with_pip09, fill = method)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = scales::label_comma()(n_genes_with_pip09)),
            vjust = -0.4, size = 2.5, na.rm = TRUE) +
  scale_fill_manual(values = method_colors) +
  labs(x = NULL, y = "Genes (PIP ≥ 0.9)", title = "Methods compared") +
  theme_masld() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 25, hjust = 1, size = 7))

opt1 <- p1a + (p1b / p1c) + plot_layout(widths = c(2, 1)) +
  plot_annotation(title = "Option 1 — Portfolio Dashboard")
save_sk(opt1, "option_01_portfolio_dashboard.pdf", w = 11, h = 6.5)

# =========================================================================
# Option 2: Causal-inference Funnel
# =========================================================================
funnel <- copy(cascade)
funnel[, step := factor(step, levels = step)]
funnel[, log_count := log10(count + 1)]
funnel[, y := -seq_len(.N)]
funnel[, half_w := log_count / max(log_count) * 4]

p2 <- ggplot(funnel) +
  geom_rect(aes(xmin = -half_w, xmax = half_w,
                ymin = y - 0.4, ymax = y + 0.4,
                fill = step), color = "white") +
  geom_text(aes(x = 0, y = y,
                label = paste0(step, " — ", scales::label_comma()(count))),
            size = 3.2, fontface = "bold") +
  scale_fill_brewer(palette = "Blues", direction = -1) +
  labs(title = "Option 2 — Causal-inference Funnel",
       subtitle = "Lead loci → finemapped → multi-ancestry → colocalized genes",
       x = NULL, y = NULL) +
  theme_void() +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, color = "gray30"),
        legend.position = "none",
        plot.margin = margin(8, 8, 8, 8))
save_sk(p2, "option_02_funnel.pdf", w = 8, h = 5)

# Option 3 (Sankey) dropped 2026-05-04: per-ancestry decomposition is not
# meaningful for cross-ancestry methods that pool EUR+EAS jointly.

# =========================================================================
# Option 4: GWAS x Method Heatmap (forensic table)
# =========================================================================
heat_long <- melt(portfolio,
                  id.vars = c("study_name","label","ancestry","trait","N_tot"),
                  measure.vars = c("susie_n_loci_converged",
                                   "susie_n_variants_pip09",
                                   "susiex_n_variants_pip09",
                                   "mesusie_n_variants_pip09",
                                   "coloc_n_genes_pp4_05"),
                  variable.name = "metric", value.name = "value")
heat_long[, metric := fcase(
  metric == "susie_n_loci_converged",  "SuSiE loci",
  metric == "susie_n_variants_pip09",  "SuSiE PIP≥0.9",
  metric == "susiex_n_variants_pip09", "SuSiE-X PIP≥0.9",
  metric == "mesusie_n_variants_pip09","meSuSiE PIP≥0.9",
  metric == "coloc_n_genes_pp4_05",    "COLOC genes")]
heat_long[, metric := factor(metric, levels = c(
  "SuSiE loci","SuSiE PIP≥0.9","SuSiE-X PIP≥0.9","meSuSiE PIP≥0.9","COLOC genes"))]
setorder(heat_long, ancestry, -N_tot)
heat_long[, label := factor(label, levels = unique(label))]

p4 <- ggplot(heat_long, aes(x = metric, y = label, fill = log10(value + 1))) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(value > 0, value, "")), size = 2.2) +
  scale_fill_gradient(low = "white", high = "#37474F", name = "log10(N+1)") +
  facet_grid(ancestry ~ ., scales = "free_y", space = "free_y") +
  labs(x = NULL, y = NULL,
       title = "Option 4 — GWAS × Method Heatmap (forensic)",
       subtitle = "Per-GWAS counts at each pipeline stage") +
  theme_masld() +
  theme(axis.text.y = element_text(size = 6),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 7),
        strip.text.y = element_text(angle = 0, face = "bold"),
        legend.position = "bottom")
save_sk(p4, "option_04_heatmap.pdf", w = 7, h = 9)

# =========================================================================
# Option 5: Concentric / Radial polar plot — retained as building block for Hybrid 12
# =========================================================================
portfolio[, x_idx := seq_len(.N)]
setorder(portfolio, ancestry, -N_tot)
portfolio[, x_idx := seq_len(.N)]

p5 <- ggplot(portfolio, aes(x = factor(x_idx), y = log10(N_tot + 1))) +
  geom_col(aes(fill = ancestry), width = 0.85, color = "white", linewidth = 0.3) +
  geom_point(aes(y = log10(N_tot + 1) + 0.3,
                 size = pmax(coloc_n_genes_pp4_05, 1)),
             color = "gray30", alpha = 0.7) +
  geom_text(aes(y = log10(N_tot + 1) - 0.4, label = label),
            size = 1.6, angle = 0, hjust = 1) +
  scale_fill_manual(values = ancestry_colors) +
  scale_size_continuous(range = c(0.5, 5), name = "COLOC\ngenes") +
  coord_polar(theta = "x", clip = "off") +
  labs(title = "Option 5 — Concentric Radial",
       subtitle = "Bars = log10(N samples), dots = #COLOC genes per GWAS") +
  theme_void() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "gray30"),
        legend.position = "right")

# =========================================================================
# Option 7 (Network) dropped 2026-05-04 — same reason as Sankey: edge weights
# to SuSiE-X / meSuSiE rely on per-GWAS attribution that doesn't exist for
# pooled cross-ancestry methods.

# =========================================================================
# Option 8: Ancestry × Trait grid (equal canvas, color + text annotation)
# Replaces area-proportional mosaic — area encoding hid AFR/SAS (small N)
# and PDFF (small UKBB MRI cohort). Numbers also re-derived to avoid
# triple-counting (Pan-UKBB ALT/AST/GGT share subjects).
# =========================================================================
cell <- portfolio[, .(n_gwas      = .N,
                      n_total     = sum(N_tot, na.rm = TRUE),
                      n_loci_susie = sum(susie_n_loci_converged, na.rm = TRUE)),
                  by = .(ancestry, trait)]
# Use UNIQUE COLOC gene counts per cell (avoids double-counting genes that hit
# multiple GWAS within the same cell). Fall back to sum if the unique table
# isn't available.
unique_genes_path <- file.path(DAT, "cell_unique_genes.csv")
if (file.exists(unique_genes_path)) {
  cg <- fread(unique_genes_path)
  cell <- merge(cell, cg, by = c("ancestry","trait"), all.x = TRUE)
  cell[, n_coloc_genes := coloc_n_genes_pp4_05_unique]
  cell[is.na(n_coloc_genes), n_coloc_genes := 0L]
  cell[, coloc_n_genes_pp4_05_unique := NULL]
} else {
  cell[, n_coloc_genes := portfolio[, sum(coloc_n_genes_pp4_05, na.rm = TRUE),
                                    by = .(ancestry, trait)][cell, on = c("ancestry","trait"), V1]]
}
cell <- cell[trait != "Other"]   # no Other-trait GWAS in registry; row dropped.
cell[, ancestry := factor(ancestry, levels = c("EUR","EAS","AFR","SAS"))]
cell[, trait := factor(trait, levels = c(
  "NAFLD/NASH","HCC","Cirrhosis","PDFF (MRI)","Liver enzymes (ALT/AST/GGT)"),
  labels = c("NAFLD/NASH","HCC","Cirrhosis","PDFF (MRI)",
             "Liver enzymes\n(ALT/AST/GGT)"))]
# Re-classify the trait factor (existing values are "Liver enzymes")
cell[, trait := fcase(
  trait == "NAFLD/NASH",                 "NAFLD/NASH",
  trait == "HCC",                        "HCC",
  trait == "Cirrhosis",                  "Cirrhosis",
  trait == "PDFF (MRI)",                 "PDFF (MRI)",
  default                                = "Liver enzymes\n(ALT/AST/GGT)")]
cell[, trait := factor(trait, levels = c(
  "NAFLD/NASH","HCC","Cirrhosis","PDFF (MRI)","Liver enzymes\n(ALT/AST/GGT)"))]

ancestry_marg <- cell[, .(n_gwas = sum(n_gwas),
                          n_total = sum(n_total),
                          n_loci = sum(n_loci_susie),
                          n_genes = sum(n_coloc_genes)),
                       by = ancestry]

# Project sunset gradient for sequential fill (matches masld_colors family).
sunset_gradient <- c("#fee292", "#fec97f", "#fcaf72", "#f9956b", "#f3796a",
                     "#e95e6d", "#db4a71", "#ca3f73", "#b83674", "#a42f74",
                     "#902972", "#7c256f")

# Project ancestry colors (masld_colors-aligned, colorblind-safe checked)
ancestry_pal <- c(EUR = "#1565C0",   # masld blue
                  EAS = "#C2185B",   # masld magenta
                  AFR = "#388E3C",   # masld green
                  SAS = "#6A1B9A")   # masld purple

p8_main <- ggplot(cell, aes(x = trait, y = ancestry)) +
  geom_tile(aes(fill = log10(n_total + 1)), color = "white", linewidth = 1.4) +
  geom_text(aes(label = paste0(n_gwas, " GWAS\nN = ",
                                scales::label_comma()(n_total), "\n",
                                n_loci_susie, " loci · ",
                                n_coloc_genes, " genes")),
            size = 3.7, lineheight = 1.0, color = "gray15") +
  scale_fill_gradientn(colors = sunset_gradient,
                       name = expression(log[10]*"(total N)"),
                       na.value = "white") +
  scale_y_discrete(limits = rev) +
  labs(x = NULL, y = NULL,
       subtitle = "Each cell: # GWAS · total sample size · finemapped loci · COLOC genes (PP4 ≥ 0.5)") +
  theme_masld() +
  theme(plot.subtitle = element_text(size = 11, color = "gray30"),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 11,
                                    lineheight = 0.9),
        axis.text.y = element_text(face = "bold", size = 13),
        legend.position = "bottom",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        legend.key.width = unit(1.2, "cm"))

p8_anc <- ggplot(ancestry_marg, aes(x = ancestry, y = n_gwas, fill = ancestry)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = paste0(n_gwas, " GWAS\n",
                                scales::label_comma()(n_total), " subj\n",
                                n_loci, " loci")),
            vjust = -0.2, size = 3.4, lineheight = 0.95) +
  scale_fill_manual(values = ancestry_pal) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.55))) +
  labs(x = NULL, y = "n GWAS", title = "Ancestry totals") +
  theme_masld() +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 13),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 11))

p8 <- (p8_main | p8_anc) + plot_layout(widths = c(2.7, 1))
save_sk(p8, "option_08_mosaic.pdf", w = 14, h = 6.5)

# =========================================================================
# Option 10: Pipeline cascade (lead loci → cross-ancestry → COLOC genes)
# Project palette: loci tier in blue family (masld_colors$control = #1565C0
# + lighter blue_gradient for sub-stages); gene tier in magenta family
# (regular ABF lighter, SuSiE COLOC darker).
# =========================================================================
cascade2 <- copy(cascade)
cascade2[, step := factor(step, levels = step)]
cascade2[, group := fcase(
  grepl("loci", step, ignore.case = TRUE), "Finemapping",
  grepl("Genes", step),                    "Colocalization")]
cascade2[, group := factor(group, levels = c("Finemapping","Colocalization"))]

# Per-step colors (within each tier, gradient by stage)
step_fill <- c(
  "Lead loci screened"                                = "#90CAF9",
  "SuSiE converged loci (canonical)"                  = "#42A5F5",
  "SuSiE-X high-PIP physical loci (cross-ancestry)"   = "#1976D2",
  "meSuSiE shared CS physical loci (cross-ancestry)"  = "#0D47A1",
  "Genes — regular (ABF) COLOC PP4 ≥ 0.5"        = "#F48FB1",
  "Genes — SuSiE COLOC PP4 ≥ 0.5"                = "#C2185B")

p10 <- ggplot(cascade2, aes(x = count, y = step, fill = step)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = scales::label_comma()(count)),
            hjust = -0.18, size = 3.8, fontface = "bold") +
  scale_x_continuous(labels = scales::label_comma(),
                     expand = expansion(mult = c(0, 0.18))) +
  scale_y_discrete(limits = rev) +
  scale_fill_manual(values = step_fill, guide = "none") +
  facet_grid(group ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(x = "Count", y = NULL,
       subtitle = "Lead loci → finemapping → cross-ancestry refinement → colocalized genes") +
  theme_masld() +
  theme(plot.subtitle = element_text(size = 11, color = "gray30"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 11),
        strip.text.y.left = element_text(angle = 0, face = "bold", size = 11),
        strip.placement = "outside",
        strip.background = element_blank(),
        panel.spacing.y = unit(0.4, "cm"))
save_sk(p10, "option_10_cascade_bar.pdf", w = 11, h = 5)

# =========================================================================
# Hybrid 11: Bubble portfolio + Cascade funnel (recommended)
# =========================================================================
portfolio[, log_N := log10(N_tot)]
portfolio[, label_short := gsub("PanUKBB ", "", label)]
hybrid_left <- ggplot(portfolio,
                      aes(x = log_N, y = reorder(label_short, log_N))) +
  geom_segment(aes(xend = 0, yend = label_short, color = ancestry), linewidth = 0.5) +
  geom_point(aes(color = ancestry, size = pmax(coloc_n_genes_pp4_05, 1))) +
  facet_grid(trait ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_size_continuous(range = c(1.2, 6), name = "COLOC genes\n(PP4 ≥ 0.5)") +
  scale_color_manual(values = ancestry_colors) +
  labs(x = "log10(N samples)", y = NULL,
       title = "GWAS portfolio") +
  theme_masld() +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 6),
        strip.text.y.left = element_text(angle = 0, face = "bold", size = 7),
        strip.placement = "outside")

hybrid_right <- p10 + labs(title = "Pipeline cascade", subtitle = NULL)

opt11 <- (hybrid_left | hybrid_right) +
  plot_layout(widths = c(1.4, 1)) +
  plot_annotation(title = "Hybrid 11 — Portfolio bubble + Cascade (RECOMMENDED)",
                  theme = theme(plot.title = element_text(face = "bold")))
save_sk(opt11, "hybrid_11_bubble_cascade.pdf", w = 13, h = 7)

# =========================================================================
# Hybrid 12: Heatmap + Radial summary inset
# =========================================================================
opt12 <- (p4 | p5) + plot_layout(widths = c(1, 1)) +
  plot_annotation(title = "Hybrid 12 — Heatmap + Radial",
                  theme = theme(plot.title = element_text(face = "bold")))
save_sk(opt12, "hybrid_12_heatmap_radial.pdf", w = 14, h = 9)

# Hybrid 13 (Bubble + Sankey) dropped 2026-05-04 along with Sankey itself.

# =========================================================================
# Hybrid 14: Cascade funnel + ancestry bar
# =========================================================================
ancestry[, ancestry := factor(ancestry, levels = c("EUR","EAS","AFR","SAS"))]
ancestry_bar <- ggplot(ancestry, aes(x = ancestry, y = n_samples_total, fill = ancestry)) +
  geom_col() +
  geom_text(aes(label = paste0(scales::label_number(suffix = "M",
                               scale = 1e-6, accuracy = 0.1)(n_samples_total),
                               "\n", n_gwas, " GWAS")),
            vjust = -0.3, size = 2.8, lineheight = 0.85) +
  scale_y_continuous(labels = scales::label_number(suffix = "M", scale = 1e-6),
                     expand = expansion(mult = c(0, 0.25))) +
  scale_fill_manual(values = ancestry_colors) +
  labs(x = NULL, y = "Cumulative samples",
       title = "Ancestry breakdown") +
  theme_masld() + theme(legend.position = "none")

ancestry_genes <- ggplot(ancestry,
                         aes(x = ancestry, y = n_coloc_genes_pp4_05, fill = ancestry)) +
  geom_col() +
  geom_text(aes(label = n_coloc_genes_pp4_05), vjust = -0.3, size = 2.8) +
  scale_fill_manual(values = ancestry_colors) +
  labs(x = NULL, y = "COLOC genes (PP4 ≥ 0.5)",
       title = "Cross-ancestry yield") +
  theme_masld() + theme(legend.position = "none") +
  expand_limits(y = max(ancestry$n_coloc_genes_pp4_05) * 1.15)

opt14 <- (p10 / (ancestry_bar | ancestry_genes)) +
  plot_layout(heights = c(1.3, 1)) +
  plot_annotation(title = "Hybrid 14 — Cascade + Ancestry bars",
                  theme = theme(plot.title = element_text(face = "bold")))
save_sk(opt14, "hybrid_14_cascade_ancestry.pdf", w = 11, h = 8)

cat("\n[render_sketches] DONE — wrote", length(list.files(OUT, pattern = "\\.pdf$")),
    "PDFs to", OUT, "\n")
