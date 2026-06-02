#!/usr/bin/env Rscript
# figS_mega_validation.R — individual panel PDFs for the mega_validation bundle.
#
# Each panel is saved as its own PDF under
#   figures/supplementary/figS_mega_validation/panels/
# with sensible standalone sizing (no composite figure — sizing tuned per panel).
#
# Panels emitted:
#   panelA_rho_heatmap.pdf            — Spearman rho on disease t-stats (4x4)
#   panelB_upset.pdf                  — UpSet of padj<0.05 DEG sets (4 methods)
#   panelC_dream_vs_eql_scatter.pdf   — dream vs edgeR-QL t-stat, mashr-colored
#   panelE_mashr_sharing.pdf          — mashr sharing-class bar
#   panelF_method_overlap_alluvial.pdf— Tier 1 method DEG overlap alluvial (4 methods)
#   panelF_method_overlap_chord.pdf   — Tier 1 method DEG overlap chord (4 methods)

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
PANEL_DIR <- file.path(PROJECT_ROOT,
  "figures/supplementary/figS_mega_validation/panels")
dir.create(PANEL_DIR, recursive = TRUE, showWarnings = FALSE)

PD <- readRDS(PANEL_RDS)

methods_pretty <- c(
  dream_t = "dream",
  vlmt    = "voomLmFit",
  dst     = "DESeq2",
  eqt     = "edgeR-QL"
)

deg_pretty <- c(
  dream = "dream",
  eq    = "edgeR-QL",
  vlm   = "voomLmFit",
  deseq = "DESeq2"
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
# Panel A — Spearman rho heatmap (4x4)
# ============================================================================
rho_dt <- as.data.table(PD$rho_t, keep.rownames = "x")
# Drop other methods before plot
drop_cols <- c("mfrt", "svt", "dlt", "elt", "trt")
rho_dt <- rho_dt[!x %in% drop_cols][, !drop_cols, with = FALSE]
rho_long <- melt(rho_dt, id.vars = "x", variable.name = "y", value.name = "rho")
rho_long[, `:=`(x = factor(methods_pretty[x], levels = methods_pretty),
                y = factor(methods_pretty[as.character(y)], levels = methods_pretty))]
pA <- ggplot(rho_long, aes(x, y, fill = rho)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.3f", rho)), size = 3) +
  scale_fill_viridis(option = "magma",
                     limits = c(min(rho_long$rho, na.rm = TRUE), 1),
                     name = expression(rho)) +
  labs(title = "Spearman rho on disease t-statistics",
       x = NULL, y = NULL) +
  panel_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(PANEL_DIR, "panelA_rho_heatmap.pdf"), pA,
       width = 5.2, height = 4.2, device = cairo_pdf)
message("Saved panelA_rho_heatmap.pdf")

# ============================================================================
# Panel B — UpSet of padj<0.05 DEG sets (4 methods)
# ============================================================================
deg_lists_raw <- PD$deg_sets[c("dream", "eq", "vlm", "deseq")]
deg_lists_raw <- deg_lists_raw[!sapply(deg_lists_raw, is.null)]
deg_lists <- deg_lists_raw[lengths(deg_lists_raw) > 0]
names(deg_lists) <- deg_pretty[names(deg_lists)]

# Hand-built UpSet matching the style of fig2_panel_nas_stage_upset.R
method_cols <- c(dream      = "#1b9e77",
                 `edgeR-QL` = "#d95f02",
                 voomLmFit  = "#7570b3",
                 DESeq2     = "#66a61e")
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
# Panel F — Method Overlap Alluvial and Chord Plots (4 methods)
# ============================================================================
# 1. Define the Tier 1 sets for each method
deg_lists_m <- list(
  dream      = PD$merged[!is.na(dream_padj) & dream_padj < 0.05 & abs(dream_logFC) > 0.5, gene],
  `edgeR-QL` = PD$merged[!is.na(eqpadj) & eqpadj < 0.05 & abs(eqlogFC) > 0.5, gene],
  voomLmFit  = PD$merged[!is.na(vlmpadj) & vlmpadj < 0.05 & abs(vlmlogFC) > 0.5, gene],
  DESeq2     = PD$merged[!is.na(dspadj) & dspadj < 0.05 & abs(dslogFC) > 0.5, gene]
)

# 2. Count how many methods each gene is called in
all_genes_m <- unique(unlist(deg_lists_m))
n_meth_per_gene <- setNames(
  rowSums(sapply(deg_lists_m, function(s) all_genes_m %in% s)),
  all_genes_m)

# 3. Create matrix of method x tier
tier_names_m <- c("Unique to method", "2 methods", "3 methods", "All 4")
mat_tier_m <- matrix(0L, length(deg_lists_m), length(tier_names_m),
                     dimnames = list(names(deg_lists_m), tier_names_m))
for (nm in names(deg_lists_m)) {
  n_others <- n_meth_per_gene[deg_lists_m[[nm]]] - 1
  for (k in 0:(length(deg_lists_m)-1)) {
    mat_tier_m[nm, k + 1] <- sum(n_others == k)
  }
}

# --- Alluvial Plot ---
df_alluv <- as.data.table(as.table(mat_tier_m))
setnames(df_alluv, c("Method", "Tier", "Freq"))
df_alluv[, Freq := as.numeric(Freq)]

method_order_top_down <- c("dream", "voomLmFit", "DESeq2", "edgeR-QL")
tier_order_top_down   <- rev(tier_names_m)

df_alluv[, Method := factor(Method, levels = method_order_top_down)]
df_alluv[, Tier   := factor(Tier,  levels = tier_order_top_down)]

method_cols <- c(dream      = "#1b9e77",
                 `edgeR-QL` = "#d95f02",
                 voomLmFit  = "#7570b3",
                 DESeq2     = "#66a61e")

tier_cols <- colorRampPalette(c("#D5D5D5", "#CE93D8", "#9575CD", "#1A237E"))(4)
names(tier_cols) <- tier_names_m

left_dt  <- df_alluv[, .(Freq = sum(Freq)), by = Method]
left_dt[, Method := factor(Method, levels = method_order_top_down)]
setorder(left_dt, Method)
total_y  <- sum(left_dt$Freq)
left_dt[, y_top := total_y - cumsum(c(0, head(Freq, -1)))]
left_dt[, y_bot := y_top - Freq]
left_dt[, y_mid := (y_top + y_bot) / 2]

right_dt <- df_alluv[, .(Freq = sum(Freq)), by = Tier]
right_dt[, Tier := factor(Tier, levels = tier_order_top_down)]
setorder(right_dt, Tier)
right_dt[, y_top := total_y - cumsum(c(0, head(Freq, -1)))]
right_dt[, y_bot := y_top - Freq]
right_dt[, y_mid := (y_top + y_bot) / 2]

p_alluv <- ggplot(df_alluv, aes(y = Freq, axis1 = Method, axis2 = Tier)) +
  geom_alluvium(aes(fill = Method), width = 1/10, alpha = 0.55,
                knot.pos = 0.42, curve_type = "sigmoid", linewidth = 0) +
  geom_stratum(aes(fill = after_stat(stratum)), width = 1/10,
               color = "white", linewidth = 0.45) +
  annotate("text", x = 1 - 0.085, y = left_dt$y_mid,
           label = as.character(left_dt$Method), hjust = 1,
           size = 2.5, fontface = "bold", color = "gray15") +
  annotate("text", x = 2 + 0.085, y = right_dt$y_mid,
           label = as.character(right_dt$Tier), hjust = 0,
           size = 2.5, fontface = "bold", color = "gray15") +
  scale_x_discrete(limits = c("Method", "Overlap"),
                   expand = c(0.30, 0.30), position = "top") +
  scale_y_continuous(expand = c(0.005, 0.005)) +
  scale_fill_manual(values = c(method_cols, tier_cols), guide = "none") +
  labs(title = "Per-method DEG overlap (Tier 1)",
       x = NULL, y = NULL) +
  theme_void() +
  theme(plot.title      = element_text(size = 9.5, face = "bold",
                                       hjust = 0.5, margin = margin(b = 6)),
        axis.text.x.top = element_text(size = 8.5, face = "bold",
                                       color = "gray15",
                                       margin = margin(b = 4)),
        plot.margin     = margin(8, 14, 6, 14))

ggsave(file.path(PANEL_DIR, "panelF_method_overlap_alluvial.pdf"), p_alluv,
       width = 6.2, height = 4.8, device = cairo_pdf)
message("Saved panelF_method_overlap_alluvial.pdf")

# --- Chord Plot ---
mat_tier_m_T <- t(mat_tier_m)
sector_cols_m <- c(tier_cols, method_cols)
gap_m <- c(setNames(rep(4, 3), tier_order_top_down[-4]),
           setNames(15, tier_order_top_down[4]),
           setNames(rep(4, 3), method_order_top_down[-4]),
           setNames(15, method_order_top_down[4]))

pdf_chord <- file.path(PANEL_DIR, "panelF_method_overlap_chord.pdf")
cairo_pdf(pdf_chord, width = 6, height = 6)
par(mar = c(1, 1, 3, 1))
circos.clear()
circos.par(start.degree = 90, canvas.xlim = c(-1.35, 1.35),
           canvas.ylim = c(-1.35, 1.35), gap.after = gap_m)
chordDiagram(
  mat_tier_m_T, grid.col = sector_cols_m, transparency = 0.35,
  annotationTrack = "grid", preAllocateTracks = list(track.height = 0.05),
  link.sort = TRUE, link.decreasing = TRUE,
  annotationTrackHeight = c(0.03, 0.03)
)

draw_label_m <- function(r_near = 1.08, r_far = 1.12, clash_deg = 22, cex = 0.70) {
  all_sec <- c(tier_order_top_down, method_order_top_down)
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
    lbl <- s
    ah <- if (th > 90 && th < 270) 1 else 0
    text(r * cos(tr) + (if (ah == 0) 0.02 else -0.02), r * sin(tr),
         lbl, cex = cex, adj = c(ah, 0.5), col = "black", xpd = TRUE)
  }
}

draw_label_m()
title(main = "Per-method DEG replication tier (Tier 1)",
      cex.main = 0.85, line = 1.5)
circos.clear()
dev.off()
message("Saved panelF_method_overlap_chord.pdf")

message("\nAll panels in: ", PANEL_DIR)
