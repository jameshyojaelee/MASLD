#!/usr/bin/env Rscript
# figS05_d_stage_da_heatmap.R
# Stage-stratified differential chromatin accessibility — bidirectional bar chart
# per cell type × contrast, showing opening (gaining) vs closing (losing) peaks.
# Hepatocytes dominate: 605 peaks in F0 vs F4, ~74% closing.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

DA_DIR <- file.path(BASE, "Analysis/ATAC/Human_Multiome/results/stage_da")
PADJ   <- 0.05

cts   <- c("Hep")
# 7 contrasts ordered: ordinal models, then fine-grained, then coarse-collapsed
tests <- c("ordinal_F0_F3_F4", "F0_vs_F3", "F3_vs_F4", "F0_vs_F4",
           "stage_ordinal",    "F0_vs_F2F3", "F2F3_vs_F4")
test_labels <- c("Ordinal F0/F3/F4",
                 "F0 vs F3",
                 "F3 vs F4",
                 "F0 vs F4",
                 "Ordinal F0/F2+F3/F4",
                 "F0 vs F2+F3",
                 "F2+F3 vs F4")

# ── load all DA results ────────────────────────────────────────────────────────
rows <- rbindlist(lapply(seq_along(tests), function(ti) {
  rbindlist(lapply(cts, function(ct) {
    f <- file.path(DA_DIR, sprintf("da_%s_%s.csv", tests[ti], ct))
    if (!file.exists(f)) return(NULL)
    d <- tryCatch(fread(f, select = c("log2FC", "padj")), error = function(e) NULL)
    if (is.null(d) || !"padj" %in% names(d)) return(NULL)
    d <- d[!is.na(padj) & padj < PADJ]
    if (nrow(d) == 0) return(NULL)
    data.table(cell_type   = ct,
               contrast    = test_labels[ti],
               n_open  = sum(d$log2FC > 0),
               n_close = sum(d$log2FC < 0))
  }), fill = TRUE)
}), fill = TRUE)

cat(sprintf("[loaded] %d cell-type × contrast combinations with DA peaks\n", nrow(rows)))

# keep only cell types with ≥ 10 total DA peaks across any contrast
ct_totals <- rows[, .(total = sum(n_open + n_close)), by = cell_type]
keep_cts  <- ct_totals[total >= 10, cell_type]
rows      <- rows[cell_type %in% keep_cts]
cat(sprintf("[filter] %d cell types with ≥10 total DA peaks: %s\n",
            length(keep_cts), paste(keep_cts, collapse = ", ")))

# melt to long for bidirectional bars
bar <- rbindlist(list(
  rows[, .(cell_type, contrast, n =  n_open,  direction = "Opening")],
  rows[, .(cell_type, contrast, n = -n_close, direction = "Closing")]
))

# factor order: contrasts left to right by biological progression
bar[, contrast   := factor(contrast, levels = test_labels)]
bar[, cell_type  := factor(cell_type, levels = keep_cts)]
bar[, direction  := factor(direction, levels = c("Opening", "Closing"))]
bar[, lbl        := fifelse(abs(n) >= 5, as.character(abs(n)), "")]

# ── plot ──────────────────────────────────────────────────────────────────────
p <- ggplot(bar, aes(x = contrast, y = n, fill = direction)) +
  geom_col(width = 0.72) +
  geom_hline(yintercept = 0, linewidth = 0.35, color = "gray60") +
  facet_wrap(~ cell_type, nrow = 1, scales = "free_y") +
  scale_fill_manual(
    values = c(Opening = masld_colors$up,    # magenta — gaining accessibility
               Closing = masld_colors$down), # blue   — losing accessibility
    name = NULL
  ) +
  scale_y_continuous(labels = function(x) comma(abs(x)),
                     expand = expansion(mult = c(0.12, 0.15))) +
  labs(x = NULL,
       y = "DA peaks  (padj < 0.05)",
       title = "Stage-stratified differential chromatin accessibility") +
  theme_masld(base_size = 10) +
  theme(strip.text       = element_text(size = 10, face = "bold"),
        axis.text.x      = element_text(size = 8, angle = 35, hjust = 1),
        axis.title.y     = element_text(size = 9),
        legend.position  = "bottom",
        legend.text      = element_text(size = 8),
        legend.key.size  = unit(0.3, "cm"),
        plot.title       = element_text(size = 11, face = "bold",
                                        margin = margin(b = 6)))

out_pdf <- file.path(FIGS05_DIR, "figS05_d_stage_da_heatmap.pdf")
ggsave(out_pdf, p,
       width  = 3 + length(keep_cts) * 1.8,
       height = 4.5,
       device = cairo_pdf)
message("Wrote: ", out_pdf)

out_csv <- sub("\\.pdf$", ".csv", out_pdf)
fwrite(bar[, .(cell_type, contrast, direction, n = abs(n))], out_csv)
message("Wrote: ", out_csv)
