#!/usr/bin/env Rscript
# figS_mega_validation_method_venn.R
# Method-overlap panels for the mega_validation bundle.
#
# Emits under figures/supplementary/figS_methods_validation/mega_validation/panels/:
#   method_venn_padj05.pdf            — 4-set Euler at padj<0.05 (dream + 3 NB/voom)
#   method_upset_tier1.pdf            — 4-set UpSet at Tier 1 (padj<0.05 & |LFC|>0.3)
#                                       style matches nas_stage_upset.pdf
#   method_venn_pairs_*.csv           — pairwise intersection summaries
#
# metafor-HKSJ removed from all comparisons (user request 2026-05-19): its
# small-K HKSJ has 0 DEGs and behaves as a heterogeneity diagnostic, not as
# a DEG voter. Atlas columns metafor_logFC / tau2 / I2 / HKSJ_padj remain
# for per-gene heterogeneity reporting.
#
# Runs in `motifbreakr` env (eulerr 7.1.0 needs RcppArmadillo compile which
# fails in rnaseq env due to SHLIB_LIBADD parsing quirk).

suppressPackageStartupMessages({
  library(data.table)
  library(eulerr)
  library(ggplot2)
  library(cowplot)
  library(scales)
  library(grid)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT  <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR <- file.path(INT, "results/mega_validation")
OUT_DIR <- file.path(BASE,
  "figures/supplementary/figS_methods_validation/mega_validation/panels")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# --- Load per-method results (metafor-HKSJ excluded from comparisons) ---
clean_id <- function(x) sub("\\..*", "", x)
load_arm <- function(path, padj_col = "padj", lfc_col = "logFC") {
  if (!file.exists(path)) return(NULL)
  dt <- fread(path)
  dt[, gene_clean := clean_id(gene)]
  list(dt = dt, padj_col = padj_col, lfc_col = lfc_col)
}

arms <- list(
  dream      = load_arm(file.path(INT, "results/integration/dream_results.csv")),
  `edgeR-QL` = load_arm(file.path(RDIR, "edgeRql/ql_results.csv")),
  voomLmFit  = load_arm(file.path(RDIR, "voomLmFit/vlm_results.csv")),
  DESeq2     = load_arm(file.path(RDIR, "DESeq2/deseq2_results.csv"))
)
arms <- arms[!sapply(arms, is.null)]
cat("Methods present:", paste(names(arms), collapse = ", "), "\n")

method_cols <- c(dream         = "#1b9e77",
                 `edgeR-QL`    = "#d95f02",
                 voomLmFit     = "#7570b3",
                 DESeq2        = "#66a61e")

# ============================================================================
# Build sets for both thresholds (all methods present)
# ============================================================================
set_padj05 <- lapply(arms, function(a) a$dt[!is.na(padj) & padj < 0.05, gene_clean])
tier1_sets <- lapply(arms, function(a)
  a$dt[!is.na(padj) & padj < 0.05 & abs(logFC) > 0.3, gene_clean])

# ============================================================================
# Panel 1 — Euler at padj < 0.05 (eulerr ellipse mode handles up to ~8 sets,
# though intersection accuracy degrades; UpSet remains the readable companion)
# ============================================================================
save_euler <- function(set_list, filename, title,
                       width = 5.5, height = 4.5) {
  fit <- euler(set_list, shape = "ellipse")
  pdf(file.path(OUT_DIR, filename), width = width, height = height)
  print(plot(fit,
             fills = list(fill = method_cols[names(set_list)], alpha = 0.55),
             edges = list(col = "white", lwd = 1.4),
             labels = list(font = 1, cex = 0.75, col = "black"),
             quantities = list(cex = 0.6, col = "black")))
  dev.off()
  cat("Saved", filename, "  (sets:",
      paste(names(set_list), lengths(set_list), sep = "=", collapse = ", "),
      ")\n")
}
message("[caption] DEG overlap at padj < 0.05")
save_euler(set_padj05, "method_venn_padj05.pdf",
           "DEG overlap at padj < 0.05")

# ============================================================================
# Panel 2 — UpSet at Tier 1 (padj<0.05 & |LFC|>0.3)
# Style follows nas_stage_upset.R: top intersection-size bar,
# bottom-right dot matrix with connecting segments, bottom-left set-size bar.
# ============================================================================
set_names <- names(tier1_sets)
N_SETS <- length(set_names)
# Cap displayed intersections — with 8 sets, 2^8-1 = 255 possible regions;
# show top by intersect size for readability.
TOP_N <- 20

# Build per-gene set membership matrix
all_genes <- unique(unlist(tier1_sets))
mem <- as.data.table(setNames(
  lapply(set_names, function(s) all_genes %in% tier1_sets[[s]]),
  paste0("in_", set_names)))
mem[, gene := all_genes]

mem_cols <- paste0("in_", set_names)
mem[, intersection_id := do.call(paste0,
       lapply(.SD, function(x) as.integer(x))),
    .SDcols = mem_cols]

intersections <- mem[, .(n_genes = .N), by = intersection_id]
intersections[, c(mem_cols) := lapply(seq_along(mem_cols), function(i) {
  substr(intersection_id, i, i) == "1"
})]
intersections[, n_sets := rowSums(.SD), .SDcols = mem_cols]
all_intersections <- intersections[n_sets > 0]
setorder(all_intersections, -n_genes)
all_intersections[, ix_rank := .I]
# Keep top-N by gene count for the figure; full table emitted to CSV
intersections <- head(all_intersections, TOP_N)
intersections[, ix_label := factor(ix_rank, levels = ix_rank)]

set_sizes <- data.table(
  method  = set_names,
  n_genes = sapply(set_names, function(s) sum(mem[[paste0("in_", s)]]))
)
set_sizes[, method := factor(method, levels = rev(set_names))]

BAR_FILL <- "#3a86ff"

# Top: intersection bar
p_top <- ggplot(intersections, aes(x = ix_label, y = n_genes)) +
  geom_col(width = 0.75, fill = BAR_FILL, color = NA) +
  geom_text(aes(label = comma(n_genes)), vjust = -0.3,
            size = 6 / ggplot2::.pt, color = "black") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.20)), labels = comma) +
  labs(x = NULL, y = "DEGs (intersection size)") +
  theme_minimal(base_size = 9) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        axis.text.x        = element_blank(),
        axis.ticks.x       = element_blank(),
        plot.margin        = margin(8, 6, 4, 6))

# Bottom-right: dot matrix with connecting segments
dot_dt <- melt(intersections[, c("ix_label", mem_cols), with = FALSE],
               id.vars = "ix_label",
               variable.name = "method", value.name = "in_set")
dot_dt[, method := factor(sub("^in_", "", method), levels = rev(set_names))]

segs <- dot_dt[in_set == TRUE,
               .(ymin = min(as.integer(method)),
                 ymax = max(as.integer(method))),
               by = ix_label]
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
  theme(panel.grid    = element_blank(),
        axis.text.x   = element_blank(),
        axis.ticks.x  = element_blank(),
        axis.text.y   = element_text(face = "plain", size = 6),
        plot.margin   = margin(0, 4, 4, 4))

# Bottom-left: set-size bar
p_left <- ggplot(set_sizes,
                 aes(y = method, x = n_genes, fill = as.character(method))) +
  geom_col(width = 0.7, color = NA) +
  geom_text(aes(label = comma(n_genes)), hjust = 1.1, size = 6 / ggplot2::.pt,
            color = "white") +
  scale_fill_manual(values = method_cols, guide = "none") +
  scale_x_reverse(expand = expansion(mult = c(0.20, 0.03)), labels = comma,
                  breaks = scales::breaks_pretty(n = 3)) +
  labs(x = "Set size", y = NULL) +
  theme_minimal(base_size = 9) +
  theme(panel.grid    = element_blank(),
        axis.text.y   = element_blank(),
        axis.ticks.y  = element_blank(),
        axis.text.x   = element_text(size = 6, angle = 35, hjust = 1),
        plot.margin   = margin(2, 4, 6, 14))

aligned <- align_plots(p_top, p_dots, align = "v", axis = "lr")
top_row    <- plot_grid(NULL, aligned[[1]], rel_widths = c(1.8, 5), nrow = 1)
bottom_row <- plot_grid(p_left, aligned[[2]], rel_widths = c(1.8, 5),
                         nrow = 1, align = "h", axis = "tb")
upset_plot <- plot_grid(top_row, bottom_row, ncol = 1,
                         rel_heights = c(1.8, 1.0))

ggsave(file.path(OUT_DIR, "method_upset_tier1.pdf"), upset_plot,
       width = 7.09, height = 3.67, device = cairo_pdf)
cat("Saved method_upset_tier1.pdf  (",
    paste(set_names, lengths(tier1_sets), sep = "=", collapse = ", "),
    ")\n")

# Full intersection table (not capped to TOP_N)
fwrite(all_intersections[, c("intersection_id", mem_cols,
                              "n_sets", "n_genes", "ix_rank"), with = FALSE],
       file.path(OUT_DIR, "method_upset_tier1_intersections.csv"))

# ============================================================================
# Pairwise summaries
# ============================================================================
pair_int <- function(sets) {
  nm <- names(sets)
  out <- expand.grid(a = nm, b = nm, stringsAsFactors = FALSE)
  out$intersect <- mapply(function(a, b) length(intersect(sets[[a]], sets[[b]])),
                          out$a, out$b)
  out$union <- mapply(function(a, b) length(union(sets[[a]], sets[[b]])),
                      out$a, out$b)
  out$jaccard <- ifelse(out$union == 0, NA_real_, out$intersect / out$union)
  setDT(out); out[]
}
fwrite(pair_int(set_padj05), file.path(OUT_DIR, "method_venn_pairs_padj05.csv"))
fwrite(pair_int(tier1_sets), file.path(OUT_DIR, "method_venn_pairs_tier1.csv"))

cat("\nDone. Panels under", OUT_DIR, "\n")
