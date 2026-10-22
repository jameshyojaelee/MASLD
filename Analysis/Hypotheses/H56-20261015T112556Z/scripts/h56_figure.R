#!/usr/bin/env Rscript
# KEY MESSAGE: for each of the 1,854 bulk calls, how far the disease log2FC moves when estimated cell mix
# enters the model, across the eight flag cells (four codings x two fraction inputs).
# Supplementary candidate panel (prespec "Destination"): x = reported unadjusted log2FC,
# y = adjusted log2FC, one vertical segment per call spanning the eight flag cells; identity and half lines.
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
OUT <- Sys.getenv("H56_OUT"); RUN <- Sys.getenv("H56_RUN")
pdf.options(useDingbats = FALSE)
gl <- fread(file.path(OUT, "call_labels_composition_sensitivity.tsv"))
gl[, lab := fifelse(is.na(label), "labels not computed",
                    fifelse(label == "at_least_half_retained_all_arms", "at least half retained, all cells",
                            fifelse(label == "more_than_half_removed_all_arms", "more than half removed, all cells",
                                    "depends on cell or not resolved")))]
cols <- c("at least half retained, all cells" = "#e14b9d", "more than half removed, all cells" = "#4baeef",
          "depends on cell or not resolved" = "#9E9E9E", "labels not computed" = "#9E9E9E")
gl[, lab := factor(lab, levels = names(cols))]
setorder(gl, lab)
cnt <- gl[, .N, by = lab]
p <- ggplot(gl) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.25, colour = "black") +
  geom_abline(slope = 0.5, intercept = 0, linewidth = 0.25, colour = "black", linetype = "dashed") +
  geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey60") +
  geom_segment(aes(x = logFC_stored, xend = logFC_stored, y = logFC_adj_min, yend = logFC_adj_max, colour = lab),
               linewidth = 0.15, alpha = 0.6) +
  scale_colour_manual(values = cols, drop = TRUE,
                      labels = function(l) paste0(l, " (n = ", cnt$N[match(l, cnt$lab)], ")")) +
  labs(x = "Reported disease log2FC, unadjusted (774 units)", y = "Disease log2FC adjusted for estimated cell mix\n(range over eight flag cells)",
       colour = NULL) +
  coord_cartesian() +
  theme_classic(base_size = 6) +
  theme(text = element_text(size = 6, colour = "black"), axis.text = element_text(size = 6, colour = "black"),
        legend.text = element_text(size = 6), legend.position = "bottom", legend.direction = "vertical",
        axis.line = element_line(linewidth = 0.25), axis.ticks = element_line(linewidth = 0.25))
dir.create(file.path(OUT, "figures"), showWarnings = FALSE)
f <- file.path(OUT, "figures", "H56_supp_unadjusted_vs_adjusted_log2FC.pdf")
ggsave(f, p, width = 3.2, height = 3.6, units = "in", useDingbats = FALSE)
message("Caption: Each vertical segment is one of the 1,854 bulk disease calls; its x position is the reported unadjusted log2FC ",
        "and it spans the adjusted log2FC across the eight flag cells (A16, C7, D1, DR-main codings x rebuilt and earlier fraction inputs). ",
        "Solid line: identity. Dashed line: half of the reported value. Adjusted values are conditional on estimated composition.")
cat("wrote", f, "\n")
