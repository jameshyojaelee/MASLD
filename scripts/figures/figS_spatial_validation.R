##############################################################################
# Supplementary Figure S14: Spatial Validation
#
# Independent spatial replication across two Visium datasets:
#   Guilliams et al. (Cell 2022) — fresh-frozen, healthy + steatotic
#   Vu et al. (JHEP Reports 2025) — FFPE, fibrosis F0–F4
#
# 7 panels:
#   (a) Violin: up-regulated DEGs are more spatially patterned
#   (b) Dot plot: which gene sets are spatially enriched
#   (c) Scatter: same genes sit in same liver zones across datasets
#   (d) Bar: periportal DEGs outnumber pericentral 2.5:1 in both
#   (e) Bar: fibrosis LR pairs dominate spatial co-expression
#   (f) Lollipop: SVG enrichment escalates with fibrosis transition
#   (g) Forest: F3->F4 pericentral enrichment (OR=8.63)
##############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(scales)
  library(ggrepel)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))

source(file.path(BASE, "scripts/figures/load_figure_data.R"))

VAL_DIR <- file.path(BASE, "Analysis/Spatial/results/validation_bulk")
PROG_DIR <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/progression")
OUT     <- file.path(FIGS05_DIR, "figS14_spatial_validation.pdf")
dir.create(FIGS05_DIR, showWarnings = FALSE, recursive = TRUE)

# Colors & labels
vu_col <- "#C2185B"
gu_col <- "#1565C0"
VU <- "Vu et al."
GU <- "Guilliams et al."
# Display labels (accession for the deposited datasets; Vu et al. has no
# accession in the registry so it is shown by author). GU/VU above remain the
# DATA KEYS that match the `dataset` column in the loaded CSVs — never relabel
# the key, only the rendered text.
GU_LAB <- "GSE192741"
VU_LAB <- VU
ds_pal <- setNames(c(gu_col, vu_col), c(GU, VU))
ds_lab <- setNames(c(GU_LAB, VU_LAB), c(GU, VU))   # display labels keyed by data value

# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------
enrich   <- fread(file.path(VAL_DIR, "deg_spatial_enrichment.csv"))
# On-disk file is conserved_core_spatial.csv (comparison=="Conserved_Core"); the
# old name conserved_spatial.csv does not exist (F199/F246). The select at line
# ~124 relabels comparison to "Conserved" so the downstream %in% filter matches.
cc_test  <- fread(file.path(VAL_DIR, "conserved_core_spatial.csv"))
merged_gu <- fread(file.path(VAL_DIR, "spatial_bulk_merged_Guilliams_et_al.csv"))
merged_vu <- fread(file.path(VAL_DIR, "spatial_bulk_merged_Vu_et_al.csv"))
zon      <- fread(file.path(VAL_DIR, "deg_zonation_combined.csv"))
lr       <- fread(file.path(VAL_DIR, "lr_pairs_combined.csv"))

# Helper: coerce Python "True"/"False" strings
fix_bool <- function(dt, cols) {
  for (col in cols) {
    if (col %in% names(dt) && is.character(dt[[col]]))
      dt[, (col) := get(col) == "True"]
  }
  dt
}

# ==========================================================================
# (a) Violin: spatial patterning of DEGs
# ==========================================================================

build_vdt <- function(dt, ds) {
  fix_bool(dt, c("is_deg", "is_up", "is_down"))
  dt[, .(I, logFC,
         group = fcase(
           is_up == TRUE,  "Up DEG",
           is_down == TRUE, "Down DEG",
           default = "Non-DEG"),
         dataset = ds)]
}

vdt <- rbind(build_vdt(merged_gu, GU), build_vdt(merged_vu, VU))
vdt[, group := factor(group, levels = c("Non-DEG", "Up DEG", "Down DEG"))]

grp_pal <- c("Non-DEG" = "#BDBDBD", "Up DEG" = vu_col, "Down DEG" = gu_col)

# P-values (vs Non-DEG, one-sided greater)
pv <- vdt[, {
  nd <- I[group == "Non-DEG"]
  up <- I[group == "Up DEG"]
  dn <- I[group == "Down DEG"]
  list(
    p_up = wilcox.test(up, nd, alternative = "greater")$p.value,
    p_dn = wilcox.test(dn, nd, alternative = "greater")$p.value
  )
}, by = dataset]

# Star labels
star <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "n.s.")))
pv[, lab_up := star(p_up)]
pv[, lab_dn := star(p_dn)]

y_top <- quantile(vdt$I, 0.995, na.rm = TRUE) * 1.2

p_a <- ggplot(vdt, aes(x = group, y = I, fill = group)) +
  geom_violin(scale = "width", width = 0.7, alpha = 0.55,
              linewidth = 0.2, color = "gray50") +
  geom_boxplot(width = 0.12, outlier.size = 0.15, outlier.alpha = 0.2,
               linewidth = 0.25, color = "gray30", fill = "white") +
  facet_wrap(~ dataset, labeller = as_labeller(ds_lab)) +
  scale_fill_manual(values = grp_pal, guide = "none") +
  geom_text(data = pv, aes(x = 2, y = y_top, label = lab_up),
            inherit.aes = FALSE, size = 2.5, color = vu_col, fontface = "bold") +
  geom_text(data = pv, aes(x = 3, y = y_top, label = lab_dn),
            inherit.aes = FALSE, size = 2.5, color = gu_col, fontface = "bold") +
  coord_cartesian(ylim = c(-0.02, y_top * 1.05)) +
  labs(x = NULL, y = "Spatial patterning\n(Moran's I)",
       title = "Up-regulated DEGs are spatially organized in liver tissue") +
  theme_masld() +
  theme(strip.text = element_text(size = 7),
        axis.text.x = element_text(size = 6, angle = 25, hjust = 1))

# ==========================================================================
# (b) Dot plot: gene set spatial enrichment
# ==========================================================================

cc_fmt <- cc_test[, .(dataset, comparison = "Conserved",
                       fold_enrichment, mannwhitney_pval)]
enrich_slim <- enrich[, .(dataset, comparison, fold_enrichment, mannwhitney_pval)]
eall <- rbind(enrich_slim, cc_fmt, fill = TRUE)
eall <- eall[comparison %in% c("Strong DEGs (|LFC|>0.5)", "Up-regulated DEGs", "Conserved")]

# Shorter labels
eall[, label := fcase(
  comparison == "Strong DEGs (|LFC|>0.5)", "Strong DEGs",
  comparison == "Up-regulated DEGs", "Up-regulated DEGs",
  comparison == "Conserved", "Conserved"
)]
eall[, label := factor(label, levels = rev(c("Strong DEGs", "Up-regulated DEGs", "Conserved")))]
eall[, sig := mannwhitney_pval < 0.001]
eall[, nlp := -log10(pmax(mannwhitney_pval, 1e-30))]

p_b <- ggplot(eall, aes(x = fold_enrichment, y = label)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray70", linewidth = 0.3) +
  geom_point(aes(color = dataset, size = nlp, shape = sig)) +
  scale_color_manual(values = ds_pal, labels = ds_lab, name = NULL) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
  scale_size_continuous(range = c(1.5, 4),
                        name = expression(-log[10] ~ italic(p)),
                        breaks = c(5, 15)) +
  scale_x_continuous(breaks = c(1, 1.5, 2, 2.5)) +
  labs(x = "Spatial enrichment (fold over non-DEGs)", y = NULL,
       title = "Spatially enriched gene sets") +
  theme_masld() +
  theme(legend.position = "right", legend.key.size = unit(0.25, "cm"))

# ==========================================================================
# (c) Scatter: cross-dataset zonation concordance
# ==========================================================================

z1 <- zon[dataset == GU, .(gene, rho1 = spearman_rho, cls1 = zonation_class)]
z2 <- zon[dataset == VU, .(gene, rho2 = spearman_rho, cls2 = zonation_class)]
zx <- merge(z1, z2, by = "gene")
zx[, agree := cls1 == cls2]

# Only highlight concordant periportal/pericentral
zx[, highlight := fcase(
  agree & cls1 == "Periportal-enriched", "Periportal",
  agree & cls1 == "Pericentral-enriched", "Pericentral",
  default = "Other"
)]

hl_pal <- c(Periportal = gu_col, Pericentral = vu_col, Other = "gray88")

# Labels: top concordant zoned genes
lab_zx <- zx[highlight != "Other"][order(-abs(rho1))][1:12]

p_c <- ggplot(zx, aes(x = rho1, y = rho2)) +
  geom_point(data = zx[highlight == "Other"], color = "gray88", size = 0.15, alpha = 0.25, shape = 16) +
  geom_point(data = zx[highlight != "Other"], aes(color = highlight),
             size = 0.8, alpha = 0.7, shape = 16) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.3) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray75", linewidth = 0.15) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray75", linewidth = 0.15) +
  geom_label_repel(data = lab_zx, aes(label = gene, color = highlight),
                   size = 1.6, max.overlaps = 20,
                   label.padding = 0.07, box.padding = 0.2,
                   segment.size = 0.1, fill = alpha("white", 0.85),
                   show.legend = FALSE) +
  scale_color_manual(values = hl_pal, name = NULL,
                     guide = guide_legend(override.aes = list(size = 2))) +
  coord_fixed() +
  labs(x = paste0("Zonation score (", GU_LAB, ")"),
       y = paste0("Zonation score (", VU_LAB, ")"),
       title = "Same genes sit in same liver zones") +
  theme_masld() +
  theme(legend.position = c(0.87, 0.15),
        legend.background = element_rect(fill = alpha("white", 0.9), color = NA),
        legend.key.size = unit(0.2, "cm"))

# ==========================================================================
# (d) Bar: periportal vs pericentral DEG counts
# ==========================================================================

zoned <- zon[zonation_class %in% c("Periportal-enriched", "Pericentral-enriched")]
zcnt <- zoned[, .N, by = .(dataset, zonation_class)]
zcnt[, zone := fifelse(grepl("Periportal", zonation_class), "PP", "PC")]
zcnt[, zone := factor(zone, levels = c("PP", "PC"))]

# Ratios
rat <- zoned[, {
  pp <- sum(zonation_class == "Periportal-enriched")
  pc <- sum(zonation_class == "Pericentral-enriched")
  list(ratio = pp / max(pc, 1))
}, by = dataset]

p_d <- ggplot(zcnt, aes(x = dataset, y = N, fill = zone)) +
  geom_col(position = position_dodge(width = 0.6), width = 0.5) +
  geom_text(aes(label = N), position = position_dodge(width = 0.6),
            vjust = -0.3, size = 2.2) +
  geom_text(data = rat,
            aes(x = dataset, y = 450,
                label = sprintf("PP:PC = %.1f:1", ratio)),
            inherit.aes = FALSE, size = 2.2, color = "gray30", fontface = "italic") +
  scale_fill_manual(values = c(PP = gu_col, PC = vu_col), name = NULL) +
  scale_x_discrete(labels = ds_lab) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "Zonated DEGs", title = "Periportal excess") +
  theme_masld() +
  theme(legend.position = c(0.85, 0.85),
        legend.background = element_rect(fill = alpha("white", 0.9), color = NA),
        axis.text.x = element_text(size = 6))

# ==========================================================================
# (e) Bar: LR spatial co-expression in Vu et al.
# ==========================================================================

lr_vu <- lr[dataset == VU]
lr_vu[, pair := paste0(ligand, "\u2013", receptor)]
lr_vu[, Category := tools::toTitleCase(category)]

setorder(lr_vu, Category, -spatial_coexpr_rho)
lr_vu[, pair := factor(pair, levels = rev(pair))]

# DEG involvement as dots
lr_vu[, deg_mark := fcase(
  ligand_is_deg & receptor_is_deg, "\u25cf\u25cf",
  ligand_is_deg | receptor_is_deg, "\u25cf",
  default = ""
)]

lr_pal <- c(Fibrosis = "#C2185B", Inflammation = "#E91E63", Sinusoidal = "#1565C0")

p_e <- ggplot(lr_vu, aes(x = spatial_coexpr_rho, y = pair, fill = Category)) +
  geom_col(width = 0.6) +
  geom_text(aes(x = -0.006, label = deg_mark), hjust = 1, size = 2, color = "gray30") +
  scale_fill_manual(values = lr_pal, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.06, 0.03))) +
  labs(x = expression("Spatial co-expression (" * rho * ")"), y = NULL,
       title = "LR pairs in MASLD tissue (Vu et al.)") +
  theme_masld() +
  theme(legend.position = c(0.82, 0.2),
        legend.background = element_rect(fill = alpha("white", 0.9), color = NA),
        legend.key.size = unit(0.25, "cm"),
        axis.text.y = element_text(size = 6))

# ==========================================================================
# (f) Lollipop: SVG enrichment escalates with fibrosis transition
# ==========================================================================

prog_svg <- tryCatch({
  cat("  Panel (f): SVG enrichment by fibrosis transition\n")
  svg_data <- fread(file.path(PROG_DIR, "spatial_validation_test1_svg_enrichment.csv"))

  # Filter to fibrosis transitions x MASLD SVGs only
  fib_svg <- svg_data[test == "transition_svg_enrichment" &
                       svg_type == "SVG_MASLD" &
                       transition %in% c("F0_to_F1", "F1_to_F2", "F2_to_F3", "F3_to_F4")]
  fib_svg[, transition := factor(
    gsub("_to_", " -> ", transition),
    levels = c("F0 -> F1", "F1 -> F2", "F2 -> F3", "F3 -> F4")
  )]
  fib_svg[, sig := padj < 0.05]
  fib_svg[, or_label := sprintf("%.2f", odds_ratio)]

  p_f <- ggplot(fib_svg, aes(x = odds_ratio, y = transition)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray60", linewidth = 0.3) +
    geom_segment(aes(x = 1, xend = odds_ratio, yend = transition,
                     color = sig), linewidth = 0.6) +
    geom_point(aes(color = sig, size = -log10(pmax(padj, 1e-30))), shape = 16) +
    geom_text(aes(label = or_label), hjust = -0.3, size = 2.3, color = "gray20") +
    scale_color_manual(values = c(`TRUE` = vu_col, `FALSE` = "gray60"),
                       labels = c(`TRUE` = "padj < 0.05", `FALSE` = "n.s."),
                       name = NULL) +
    scale_size_continuous(range = c(2, 4.5),
                          name = expression(-log[10] ~ padj),
                          breaks = c(2, 7, 13)) +
    scale_x_continuous(limits = c(0.8, 3.3)) +
    labs(x = "Odds ratio (SVG enrichment)", y = NULL,
         title = "SVG enrichment escalates with fibrosis stage") +
    theme_masld() +
    theme(legend.position = "right",
          legend.key.size = unit(0.25, "cm"))

  p_f
}, error = function(e) {
  cat("    Panel (f) FAILED:", conditionMessage(e), "\n")
  placeholder("Panel (f) failed")
})

# ==========================================================================
# (g) Forest: pericentral enrichment by fibrosis transition
# ==========================================================================

prog_zon <- tryCatch({
  cat("  Panel (g): Pericentral enrichment by fibrosis transition\n")
  zon_data <- fread(file.path(PROG_DIR, "spatial_validation_test2_zonation.csv"))

  # Filter to transition zonation enrichment for fibrosis + pericentral
  fib_zon <- zon_data[test == "transition_zonation_enrichment" &
                       zone == "Pericentral" &
                       gene_set %in% c("F0_to_F1", "F1_to_F2", "F2_to_F3", "F3_to_F4")]
  fib_zon[, transition := factor(
    gsub("_to_", " -> ", gene_set),
    levels = c("F0 -> F1", "F1 -> F2", "F2 -> F3", "F3 -> F4")
  )]
  fib_zon[, sig := padj < 0.05]
  fib_zon[, or_label := sprintf("%.1f", odds_ratio)]

  # Also add S1 pericentral enrichment
  s1_zon <- zon_data[test == "zonation_enrichment" &
                      zone == "Pericentral" &
                      gene_set == "divergence_S1_up"]
  if (nrow(s1_zon) > 0) {
    s1_zon[, transition := "S1 subtype"]
    s1_zon[, sig := padj < 0.05]
    s1_zon[, or_label := sprintf("%.1f", odds_ratio)]
    fib_zon <- rbind(fib_zon, s1_zon, fill = TRUE)
    fib_zon[, transition := factor(transition,
      levels = c("F0 -> F1", "F1 -> F2", "F2 -> F3", "F3 -> F4", "S1 subtype"))]
  }

  p_g <- ggplot(fib_zon, aes(x = odds_ratio, y = transition)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray60", linewidth = 0.3) +
    geom_segment(aes(x = 1, xend = odds_ratio, yend = transition,
                     color = sig), linewidth = 0.6) +
    geom_point(aes(color = sig, size = -log10(pmax(padj, 1e-30))), shape = 16) +
    geom_text(aes(label = or_label), hjust = -0.3, size = 2.3, color = "gray20") +
    scale_color_manual(values = c(`TRUE` = vu_col, `FALSE` = "gray60"),
                       labels = c(`TRUE` = "padj < 0.05", `FALSE` = "n.s."),
                       name = NULL) +
    scale_size_continuous(range = c(2, 4.5),
                          name = expression(-log[10] ~ padj),
                          breaks = c(2, 5, 12)) +
    scale_x_continuous(limits = c(0.5, 10.5)) +
    labs(x = "Odds ratio (pericentral enrichment)", y = NULL,
         title = "F3->F4 genes concentrate in pericentral zone") +
    theme_masld() +
    theme(legend.position = "right",
          legend.key.size = unit(0.25, "cm"))

  p_g
}, error = function(e) {
  cat("    Panel (g) FAILED:", conditionMessage(e), "\n")
  placeholder("Panel (g) failed")
})

# ==========================================================================
# Compose
# ==========================================================================
cat("Composing figure...\n")

fig <- (p_a | p_b) /
       (p_c | p_d) /
       (p_e) /
       (prog_svg | prog_zon) +
  plot_layout(heights = c(1, 1.2, 0.8, 0.8)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 9, face = "bold"))

save_fig_tall(fig, OUT, width = fig_full_width, height = 10.5)
cat("Saved:", OUT, "\n")

# Final panels — individual panel PDFs
PANELS05 <- FIGS05_DIR
save_fig(p_a, file.path(PANELS05, "spatial_val_a_violin.pdf"), width = fig_full_width, height = 2.5)
save_fig(p_b, file.path(PANELS05, "spatial_val_b_enrichment.pdf"), width = fig_half_width, height = 2)
save_fig(p_c, file.path(PANELS05, "spatial_val_c_zonation.pdf"), width = fig_half_width, height = 3.5)
save_fig(p_d, file.path(PANELS05, "spatial_val_d_pp_pc.pdf"), width = fig_half_width, height = 2.5)
save_fig(p_e, file.path(PANELS05, "spatial_val_e_lr.pdf"), width = fig_full_width, height = 2.5)
save_fig(prog_svg, file.path(PANELS05, "spatial_val_f_svg_escalation.pdf"), width = fig_half_width, height = 2.5)
save_fig(prog_zon, file.path(PANELS05, "spatial_val_g_pericentral.pdf"), width = fig_half_width, height = 2.5)

cat("Done.\n")
