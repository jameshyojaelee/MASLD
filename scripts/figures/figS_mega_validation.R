#!/usr/bin/env Rscript
# figS_mega_validation.R — individual panel PDFs for the mega_validation bundle.
#
# Each panel is saved as its own PDF under
#   figures/supplementary/figS_methods_validation/mega_validation/panels/
# with sensible standalone sizing (no composite figure — sizing tuned per panel).
#
# Panels emitted:
#   panelA_rho_heatmap.pdf            — Spearman rho on logFC (6x6, 2 metafor engines)
#   panelB_upset.pdf                  — UpSet of padj<0.05 DEG sets (6 methods)
#   panelC_dream_vs_eql_scatter.pdf   — dream vs edgeR-QL t-stat, mashr-colored
#   panelC2_dream_vs_meta_scatter.pdf — dream vs metafor(voom) logFC, I2-colored
#   panelE_mashr_sharing.pdf          — mashr sharing-class bar
#   panelF_method_overlap_alluvial.pdf— Tier 1 method DEG overlap alluvial (6 methods)
#   panelF_method_overlap_chord.pdf   — Tier 1 method DEG overlap chord (6 methods)

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(viridis)
  library(cowplot); library(scales); library(grid)
  library(ggalluvial); library(circlize)
})

PROJECT_ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT <- file.path(PROJECT_ROOT, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
PANEL_RDS <- file.path(INT,
  "results/mega_validation/concordance/panel_data_A_to_F.rds")
META_ENGINES <- file.path(INT,
  "results/integration/multimethod_validation/sensitivity/deseq2_vs_limma_metafor_comparison.csv")
PANEL_DIR <- file.path(PROJECT_ROOT,
  "figures/supplementary/figS_methods_validation/mega_validation/panels")
dir.create(PANEL_DIR, recursive = TRUE, showWarnings = FALSE)

PD <- readRDS(PANEL_RDS)

# ── Load TWO metafor engine arms and merge into PD$merged ──────────────
# Sensitivity run (binary estimand) computes metafor on BOTH engines identically:
#   metafor(voom)   = limma-voom -> metafor   (lv_logFC / lv_padj)
#   metafor(DESeq2) = DESeq2     -> metafor   (ds_logFC / ds_padj)
# Gene IDs are version-suffixed and match PD$merged$gene DIRECTLY (no stripping).
meta_engines <- fread(META_ENGINES)
merged <- copy(PD$merged)
merged <- merge(merged,
                meta_engines[, .(gene, lv_logFC, lv_padj, ds_logFC, ds_padj)],
                by = "gene", all.x = TRUE)
cat(sprintf("metafor(voom)   merged: %d / %d genes matched\n",
            sum(!is.na(merged$lv_logFC)), nrow(merged)))
cat(sprintf("metafor(DESeq2) merged: %d / %d genes matched\n",
            sum(!is.na(merged$ds_logFC)), nrow(merged)))

methods_pretty <- c(
  dream_t  = "dream",
  vlmt     = "voomLmFit",
  dst      = "DESeq2",
  eqt      = "edgeR-QL",
  lv_lfc   = "metafor(voom)",
  ds_lfc   = "metafor(DESeq2)"
)

deg_pretty <- c(
  dream    = "dream",
  eq       = "edgeR-QL",
  vlm      = "voomLmFit",
  deseq    = "DESeq2",
  metaLV   = "metafor(voom)",
  metaDS   = "metafor(DESeq2)"
)

method_cols <- c(
  dream              = "#1b9e77",
  `edgeR-QL`         = "#d95f02",
  voomLmFit          = "#7570b3",
  DESeq2             = "#66a61e",
  `metafor(voom)`    = "#1f78b4",
  `metafor(DESeq2)`  = "#a6cee3"
)

share_cols <- c(pan_cohort = "#1b9e77",
                majority_shared = "#66a61e",
                cohort_specific = "#e7298a",
                divergent = "#e6ab02",
                null = "#9E9E9E")

panel_theme <- theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11),
        plot.margin = margin(8, 10, 8, 10))

# ============================================================================
# Panel A — Spearman rho heatmap (6x6, 2 metafor engines)
# ============================================================================
# Compute rho on logFC (same scale for all 6 methods; t-stats not available
# for the metafor engines so we use logFC throughout for consistency).
lfc_cols <- c(dream = "dream_logFC", `edgeR-QL` = "eqlogFC",
              voomLmFit = "vlmlogFC", DESeq2 = "dslogFC",
              `metafor(voom)` = "lv_logFC", `metafor(DESeq2)` = "ds_logFC")
lfc_mat <- as.matrix(merged[, lfc_cols, with = FALSE])
colnames(lfc_mat) <- names(lfc_cols)
rho_5 <- cor(lfc_mat, method = "spearman", use = "pairwise.complete.obs")
cat(sprintf("rho(metafor(voom), metafor(DESeq2)) = %.4f\n",
            rho_5["metafor(voom)", "metafor(DESeq2)"]))

rho_long5 <- as.data.table(rho_5, keep.rownames = "x")
rho_long5 <- melt(rho_long5, id.vars = "x", variable.name = "y", value.name = "rho")
method_order5 <- names(lfc_cols)
rho_long5[, x := factor(x, levels = method_order5)]
rho_long5[, y := factor(y, levels = method_order5)]

pA <- ggplot(rho_long5, aes(x, y, fill = rho)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.3f", rho)), size = 2.8) +
  scale_fill_viridis(option = "magma",
                     limits = c(min(rho_long5$rho, na.rm = TRUE), 1),
                     name = "Spearman rho") +
  labs(title = "Spearman rho on logFC (6 methods, 2 metafor engines)",
       x = NULL, y = NULL) +
  panel_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(PANEL_DIR, "panelA_rho_heatmap.pdf"), pA,
       width = 6.0, height = 5.0, device = cairo_pdf)
message("Saved panelA_rho_heatmap.pdf")

# ============================================================================
# Panel B — UpSet of padj<0.05 DEG sets (6 methods, 2 metafor engines)
# ============================================================================
metaLV_degs <- merged[!is.na(lv_padj) & lv_padj < 0.05, gene]
metaDS_degs <- merged[!is.na(ds_padj) & ds_padj < 0.05, gene]
deg_lists_raw <- c(PD$deg_sets[c("dream", "eq", "vlm", "deseq")],
                   list(metaLV = metaLV_degs, metaDS = metaDS_degs))
deg_lists_raw <- deg_lists_raw[!sapply(deg_lists_raw, is.null)]
deg_lists <- deg_lists_raw[lengths(deg_lists_raw) > 0]
names(deg_lists) <- deg_pretty[names(deg_lists)]
set_names <- names(deg_lists)
all_genes <- unique(unlist(deg_lists))
mem <- as.data.table(setNames(
  lapply(set_names, function(s) all_genes %in% deg_lists[[s]]),
  paste0("in_", set_names)))
mem[, gene := all_genes]
mem_cols <- paste0("in_", set_names)
mem[, intersection_id := do.call(paste0,
       lapply(.SD, function(x) as.integer(x))), .SDcols = mem_cols]

intersections <- mem[, .(n_genes = .N), by = intersection_id]
intersections[, c(mem_cols) := lapply(seq_along(mem_cols), function(i) {
  substr(intersection_id, i, i) == "1"
})]
intersections[, n_sets := rowSums(.SD), .SDcols = mem_cols]
all_intersections <- intersections[n_sets > 0]
setorder(all_intersections, -n_genes)
all_intersections[, ix_rank := .I]

intersections <- all_intersections
intersections[, ix_label := factor(ix_rank, levels = ix_rank)]

set_sizes <- data.table(
  method  = set_names,
  n_genes = sapply(set_names, function(s) sum(mem[[paste0("in_", s)]]))
)
set_sizes[, method := factor(method, levels = rev(set_names))]

BAR_FILL <- "#3a86ff"
p_top <- ggplot(intersections, aes(x = ix_label, y = n_genes)) +
  geom_col(width = 0.75, fill = BAR_FILL, color = NA) +
  geom_text(aes(label = comma(n_genes)), vjust = -0.3,
            size = 2.2, color = "gray20") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.20)), labels = comma) +
  labs(x = NULL, y = "DEGs (intersection size)") +
  theme_minimal(base_size = 9) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        axis.text.x        = element_blank(),
        axis.ticks.x       = element_blank(),
        plot.margin        = margin(8, 6, 4, 6))

dot_dt <- melt(intersections[, c("ix_label", mem_cols), with = FALSE],
               id.vars = "ix_label",
               variable.name = "method", value.name = "in_set")
dot_dt[, method := factor(sub("^in_", "", method), levels = rev(set_names))]
segs <- dot_dt[in_set == TRUE,
               .(ymin = min(as.integer(method)),
                 ymax = max(as.integer(method))), by = ix_label]
segs <- segs[ymin != ymax]

p_dots <- ggplot(dot_dt, aes(x = ix_label, y = method)) +
  geom_point(aes(color = in_set), size = 2.2) +
  { if (nrow(segs)) geom_segment(data = segs,
      aes(x = ix_label, xend = ix_label, y = ymin, yend = ymax),
      inherit.aes = FALSE, color = "#212121", linewidth = 0.5) } +
  scale_color_manual(values = c("TRUE" = "#212121", "FALSE" = "#E0E0E0"),
                     guide = "none") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(),
        axis.text.x   = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y   = element_text(face = "bold", size = 8),
        plot.margin   = margin(0, 4, 4, 4))

p_left <- ggplot(set_sizes,
                 aes(y = method, x = n_genes, fill = as.character(method))) +
  geom_col(width = 0.7, color = NA) +
  geom_text(aes(label = comma(n_genes)), hjust = 1.1, size = 2.2,
            color = "white") +
  scale_fill_manual(values = method_cols, guide = "none") +
  scale_x_reverse(expand = expansion(mult = c(0.20, 0.03)), labels = comma,
                  breaks = scales::breaks_pretty(n = 3)) +
  labs(x = "Set size", y = NULL) +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(),
        axis.text.y  = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x  = element_text(size = 7, angle = 35, hjust = 1),
        plot.margin  = margin(2, 4, 6, 14))

aligned <- align_plots(p_top, p_dots, align = "v", axis = "lr")
top_row    <- plot_grid(NULL, aligned[[1]], rel_widths = c(1.8, 5), nrow = 1)
bottom_row <- plot_grid(p_left, aligned[[2]], rel_widths = c(1.8, 5),
                         nrow = 1, align = "h", axis = "tb")
upset_padj05 <- plot_grid(top_row, bottom_row, ncol = 1,
                           rel_heights = c(1.8, 1.0))
ggsave(file.path(PANEL_DIR, "panelB_upset.pdf"), upset_padj05,
       width = 8.5, height = 4.4, device = cairo_pdf)
fwrite(all_intersections[, c("intersection_id", mem_cols, "n_sets",
                          "n_genes", "ix_rank"), with = FALSE],
       file.path(PANEL_DIR, "panelB_upset_intersections.csv"))
message("Saved panelB_upset.pdf (hand-built, padj<0.05)")

# ============================================================================
# Panel C — dream vs edgeR-QL t-stat scatter (mashr-colored)
# ============================================================================
sc <- PD$merged
mashr_has_data <- "mashr_sharing_class" %in% names(sc) &&
                  any(!is.na(sc$mashr_sharing_class))

pC <- if (mashr_has_data) {
  ggplot(sc[!is.na(eqt)],
         aes(dream_t, eqt, color = mashr_sharing_class)) +
    geom_point(size = 0.6, alpha = 0.4) +
    geom_abline(linetype = "dashed", color = "#444444") +
    scale_color_manual(values = share_cols, name = "mashr class",
                       na.value = "#cccccc") +
    labs(title = "dream vs edgeR-QL t-statistic",
         x = "dream t", y = "edgeR-QL signed-sqrt(F)") +
    guides(color = guide_legend(override.aes = list(size = 2.4, alpha = 1))) +
    panel_theme
} else {
  ggplot(sc[!is.na(eqt)], aes(dream_t, eqt)) +
    geom_point(size = 0.6, alpha = 0.4, color = "#3a86ff") +
    geom_abline(linetype = "dashed", color = "#444444") +
    labs(title = "dream vs edgeR-QL t-statistic",
         x = "dream t", y = "edgeR-QL signed-sqrt(F)") +
    panel_theme
}
ggsave(file.path(PANEL_DIR, "panelC_dream_vs_eql_scatter.pdf"), pC,
       width = 6, height = 5.5, device = cairo_pdf)
message("Saved panelC_dream_vs_eql_scatter.pdf")

# ============================================================================
# Panel C2 — dream vs metafor(voom) logFC scatter (I2-coloured)
# ============================================================================
lim_c2 <- max(abs(c(merged$dream_logFC, merged$lv_logFC)), na.rm = TRUE) * 1.05
rho_c2 <- cor(merged$dream_logFC, merged$lv_logFC,
              method = "spearman", use = "complete.obs")
pC2 <- ggplot(merged[!is.na(lv_logFC)],
              aes(dream_logFC, lv_logFC, colour = I2)) +
  geom_point(size = 0.3, alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "#B0BEC5", linewidth = 0.4) +
  scale_colour_viridis_c(name = "I2 (%)", limits = c(0, 100),
                         option = "plasma", direction = -1,
                         na.value = "#CFD8DC") +
  coord_fixed(xlim = c(-lim_c2, lim_c2), ylim = c(-lim_c2, lim_c2)) +
  labs(title = sprintf("dream vs metafor(voom) logFC  |  rho = %.3f", rho_c2),
       x = "dream logFC", y = "metafor(voom) logFC") +
  panel_theme +
  theme(legend.position = "right")
ggsave(file.path(PANEL_DIR, "panelC2_dream_vs_meta_scatter.pdf"), pC2,
       width = 6, height = 5.5, device = cairo_pdf)
message("Saved panelC2_dream_vs_meta_scatter.pdf")

# ============================================================================
# Panel E — mashr sharing-class bar
# ============================================================================
if (nrow(PD$sharing_classes) > 0) {
  sc_classes <- PD$sharing_classes[, .N, by = sharing_class][order(-N)]
  sc_classes[, pct := 100 * N / sum(N)]
  pE <- ggplot(sc_classes,
               aes(reorder(sharing_class, N), N, fill = sharing_class)) +
    geom_col() + coord_flip() +
    geom_text(aes(label = sprintf("%d (%.1f%%)", N, pct)),
              hjust = -0.05, size = 3.2) +
    scale_fill_manual(values = share_cols, guide = "none") +
    labs(title = "mashr sharing-class distribution",
         x = NULL, y = "Genes") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
    panel_theme
  ggsave(file.path(PANEL_DIR, "panelE_mashr_sharing.pdf"), pE,
         width = 6, height = 3.5, device = cairo_pdf)
  message("Saved panelE_mashr_sharing.pdf")
} else {
  message("Skipped panelE_mashr_sharing.pdf — no mashr data")
}

# ============================================================================
# Panel F — Method Overlap Alluvial (6 methods, 2 metafor engines); two thresholds
# ============================================================================
method_order_top_down <- c("dream", "voomLmFit", "DESeq2", "edgeR-QL",
                           "metafor(voom)", "metafor(DESeq2)")

make_alluvial <- function(deg_lists, title, subtitle,
                          outfile, width = 6.2, height = 4.8) {
  n_m        <- length(deg_lists)
  tier_names <- c("Unique to method",
                  paste(2:(n_m - 1), "methods"),
                  sprintf("All %d", n_m))

  all_g <- unique(unlist(deg_lists))
  n_per <- setNames(rowSums(sapply(deg_lists, function(s) all_g %in% s)), all_g)

  mat <- matrix(0L, n_m, n_m,
                dimnames = list(names(deg_lists), tier_names))
  for (nm in names(deg_lists)) {
    others <- n_per[deg_lists[[nm]]] - 1
    for (k in 0:(n_m - 1)) mat[nm, k + 1] <- sum(others == k)
  }

  tc <- colorRampPalette(c("#D5D5D5", "#CE93D8", "#9575CD", "#4527A0", "#1A237E"))(n_m)
  names(tc) <- tier_names
  tier_ord <- rev(tier_names)

  df <- as.data.table(as.table(mat))
  setnames(df, c("Method", "Tier", "Freq"))
  df[, Freq   := as.numeric(Freq)]
  df[, Method := factor(Method, levels = method_order_top_down)]
  df[, Tier   := factor(Tier,   levels = tier_ord)]

  ld <- df[, .(Freq = sum(Freq)), by = Method]
  ld[, Method := factor(Method, levels = method_order_top_down)]
  setorder(ld, Method)
  tot <- sum(ld$Freq)
  ld[, y_top := tot - cumsum(c(0, head(Freq, -1)))]
  ld[, y_mid := y_top - Freq / 2]

  rd <- df[, .(Freq = sum(Freq)), by = Tier]
  rd[, Tier := factor(Tier, levels = tier_ord)]
  setorder(rd, Tier)
  rd[, y_top := tot - cumsum(c(0, head(Freq, -1)))]
  rd[, y_mid := y_top - Freq / 2]

  p <- ggplot(df, aes(y = Freq, axis1 = Method, axis2 = Tier)) +
    geom_alluvium(aes(fill = Method), width = 1/10, alpha = 0.55,
                  knot.pos = 0.42, curve_type = "sigmoid", linewidth = 0) +
    geom_stratum(aes(fill = after_stat(stratum)), width = 1/10,
                 color = "white", linewidth = 0.45) +
    annotate("text", x = 1 - 0.085, y = ld$y_mid,
             label = as.character(ld$Method), hjust = 1,
             size = 2.5, fontface = "bold", color = "gray15") +
    annotate("text", x = 2 + 0.085, y = rd$y_mid,
             label = as.character(rd$Tier), hjust = 0,
             size = 2.5, fontface = "bold", color = "gray15") +
    scale_x_discrete(limits = c("Method", "Overlap"),
                     expand = c(0.30, 0.30), position = "top") +
    scale_y_continuous(expand = c(0.005, 0.005)) +
    scale_fill_manual(values = c(method_cols, tc), guide = "none") +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_void() +
    theme(plot.title    = element_text(size = 9.5, face = "bold",
                                       hjust = 0.5, margin = margin(b = 2)),
          plot.subtitle = element_text(size = 8, hjust = 0.5,
                                       color = "gray40", margin = margin(b = 6)),
          axis.text.x.top = element_text(size = 8.5, face = "bold",
                                         color = "gray15", margin = margin(b = 4)),
          plot.margin   = margin(8, 14, 6, 14))

  ggsave(file.path(PANEL_DIR, outfile), p,
         width = width, height = height, device = cairo_pdf)
  message("Saved ", outfile)
}

# panelF — Tier 1: padj < 0.05 AND |logFC| > 0.5
make_alluvial(
  deg_lists = list(
    dream             = merged[!is.na(dream_padj) & dream_padj < 0.05 & abs(dream_logFC) > 0.5, gene],
    `edgeR-QL`        = merged[!is.na(eqpadj)    & eqpadj    < 0.05 & abs(eqlogFC)    > 0.5, gene],
    voomLmFit         = merged[!is.na(vlmpadj)   & vlmpadj   < 0.05 & abs(vlmlogFC)   > 0.5, gene],
    DESeq2            = merged[!is.na(dspadj)    & dspadj    < 0.05 & abs(dslogFC)    > 0.5, gene],
    `metafor(voom)`   = merged[!is.na(lv_padj)   & lv_padj   < 0.05 & abs(lv_logFC)   > 0.5, gene],
    `metafor(DESeq2)` = merged[!is.na(ds_padj)   & ds_padj   < 0.05 & abs(ds_logFC)   > 0.5, gene]
  ),
  title    = "Per-method DEG overlap (Tier 1)",
  subtitle = "padj < 0.05  |  |logFC| > 0.5",
  outfile  = "panelF_method_overlap_alluvial.pdf"
)

# panelF2 — padj < 0.05 only (no LFC cutoff)
make_alluvial(
  deg_lists = list(
    dream             = merged[!is.na(dream_padj) & dream_padj < 0.05, gene],
    `edgeR-QL`        = merged[!is.na(eqpadj)    & eqpadj    < 0.05, gene],
    voomLmFit         = merged[!is.na(vlmpadj)   & vlmpadj   < 0.05, gene],
    DESeq2            = merged[!is.na(dspadj)    & dspadj    < 0.05, gene],
    `metafor(voom)`   = merged[!is.na(lv_padj)   & lv_padj   < 0.05, gene],
    `metafor(DESeq2)` = merged[!is.na(ds_padj)   & ds_padj   < 0.05, gene]
  ),
  title    = "Per-method DEG overlap",
  subtitle = "padj < 0.05  (no LFC cutoff)",
  outfile  = "panelF2_method_overlap_alluvial_padj05.pdf"
)

# --- Chord Plot (Tier 1 only) ---
deg_lists_chord <- list(
  dream             = merged[!is.na(dream_padj) & dream_padj < 0.05 & abs(dream_logFC) > 0.5, gene],
  `edgeR-QL`        = merged[!is.na(eqpadj)    & eqpadj    < 0.05 & abs(eqlogFC)    > 0.5, gene],
  voomLmFit         = merged[!is.na(vlmpadj)   & vlmpadj   < 0.05 & abs(vlmlogFC)   > 0.5, gene],
  DESeq2            = merged[!is.na(dspadj)    & dspadj    < 0.05 & abs(dslogFC)    > 0.5, gene],
  `metafor(voom)`   = merged[!is.na(lv_padj)   & lv_padj   < 0.05 & abs(lv_logFC)   > 0.5, gene],
  `metafor(DESeq2)` = merged[!is.na(ds_padj)   & ds_padj   < 0.05 & abs(ds_logFC)   > 0.5, gene]
)
n_m         <- length(deg_lists_chord)
tier_names_chord <- c("Unique to method", paste(2:(n_m-1), "methods"),
                      sprintf("All %d", n_m))
all_g_c     <- unique(unlist(deg_lists_chord))
n_per_c     <- setNames(rowSums(sapply(deg_lists_chord, function(s) all_g_c %in% s)), all_g_c)
mat_chord   <- matrix(0L, n_m, n_m,
                      dimnames = list(names(deg_lists_chord), tier_names_chord))
for (nm in names(deg_lists_chord)) {
  others <- n_per_c[deg_lists_chord[[nm]]] - 1
  for (k in 0:(n_m-1)) mat_chord[nm, k+1] <- sum(others == k)
}
tier_ord_chord <- rev(tier_names_chord)
tc_chord <- colorRampPalette(c("#D5D5D5", "#CE93D8", "#9575CD", "#4527A0", "#1A237E"))(n_m)
names(tc_chord) <- tier_names_chord
sector_cols_m <- c(tc_chord, method_cols)
gap_m <- c(setNames(rep(4, n_m-1), tier_ord_chord[-n_m]),
           setNames(15, tier_ord_chord[n_m]),
           setNames(rep(4, n_m-1), method_order_top_down[-n_m]),
           setNames(15, method_order_top_down[n_m]))

pdf_chord <- file.path(PANEL_DIR, "panelF_method_overlap_chord.pdf")
cairo_pdf(pdf_chord, width = 6, height = 6)
par(mar = c(1, 1, 3, 1))
circos.clear()
circos.par(start.degree = 90, canvas.xlim = c(-1.35, 1.35),
           canvas.ylim = c(-1.35, 1.35), gap.after = gap_m)
chordDiagram(
  t(mat_chord), grid.col = sector_cols_m, transparency = 0.35,
  annotationTrack = "grid", preAllocateTracks = list(track.height = 0.05),
  link.sort = TRUE, link.decreasing = TRUE,
  annotationTrackHeight = c(0.03, 0.03)
)
draw_label_m <- function(r_near = 1.08, r_far = 1.12, clash_deg = 22, cex = 0.70) {
  all_sec <- c(tier_ord_chord, method_order_top_down)
  si <- lapply(all_sec, function(s) {
    set.current.cell(sector.index = s, track.index = 1)
    th <- circlize(mean(get.cell.meta.data("xlim")), 1,
                   sector.index = s, track.index = 1)[1, "theta"]
    list(s = s, th = (th + 360) %% 360)
  })
  si <- si[order(sapply(si, `[[`, "th"), decreasing = TRUE)]
  rails <- rep("near", length(si))
  for (i in seq_along(si)[-1]) {
    d <- min(abs(si[[i]]$th - si[[i-1]]$th), 360 - abs(si[[i]]$th - si[[i-1]]$th))
    if (d < clash_deg) rails[i] <- if (rails[i-1] == "near") "far" else "near"
  }
  for (i in seq_along(si)) {
    s <- si[[i]]$s; th <- si[[i]]$th
    r  <- if (rails[i] == "near") r_near else r_far
    tr <- th * pi / 180
    ah <- if (th > 90 && th < 270) 1 else 0
    text(r * cos(tr) + (if (ah == 0) 0.02 else -0.02), r * sin(tr),
         s, cex = cex, adj = c(ah, 0.5), col = "black", xpd = TRUE)
  }
}
draw_label_m()
title(main = "Per-method DEG replication tier (padj<0.05, |logFC|>0.5)",
      cex.main = 0.80, line = 1.5)
circos.clear()
dev.off()
message("Saved panelF_method_overlap_chord.pdf")

# ============================================================================
# Panel H — Robustness bar: dream Tier-1 DEGs confirmed by other methods
# ============================================================================
ROB_CSV <- file.path(INT, "results/mega_validation/concordance/robustness_summary.csv")
rob <- fread(ROB_CSV)

# Human-readable labels and a logical order (best → worst confirmation)
flag_labels <- c(
  confirmed_all      = "Confirmed by all methods",
  confirmed_majority = "Confirmed by majority",
  confirmed_minority = "Confirmed by minority",
  primary_only       = "dream only (no confirmation)",
  discordant         = "Discordant"
)
flag_cols <- c(
  confirmed_all      = "#1b9e77",
  confirmed_majority = "#66a61e",
  confirmed_minority = "#e6ab02",
  primary_only       = "#d95f02",
  discordant         = "#e31a1c"
)
flag_order <- names(flag_labels)

rob[, label := flag_labels[flag]]
rob[, label := factor(label, levels = rev(flag_labels[flag_order]))]
rob[, fill  := flag_cols[flag]]
rob[, pct_label := sprintf("%d  (%.1f%%)", n, pct_tier1)]

n_tier1 <- sum(rob$n)

pH <- ggplot(rob, aes(x = n, y = label, fill = flag)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = pct_label), hjust = -0.08, size = 3.0, color = "gray20") +
  scale_fill_manual(values = flag_cols, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.22)),
                     labels = scales::comma) +
  labs(
    title    = "dream Tier-1 DEG confirmation across methods",
    subtitle = sprintf("n = %s dream Tier-1 DEGs (padj < 0.05, |logFC| > 0.5)  |  voters: edgeR-QL, voomLmFit, DESeq2",
                       scales::comma(n_tier1)),
    x = "Number of genes", y = NULL
  ) +
  panel_theme +
  theme(axis.text.y = element_text(size = 9.5))

ggsave(file.path(PANEL_DIR, "panelH_robustness_bar.pdf"), pH,
       width = 7.5, height = 3.2, device = cairo_pdf)
message("Saved panelH_robustness_bar.pdf")

message("\nAll panels in: ", PANEL_DIR)
