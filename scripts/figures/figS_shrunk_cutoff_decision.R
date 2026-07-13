#!/usr/bin/env Rscript
# Decision panel: known-MASLD recall + DEG set-size vs ashr-shrunk |log2FC| cutoff.
# Single instrument for choosing the canonical shrunk LFC floor.
# Shrunk gate = lfsr < 0.05 (the canonical DEG definition); cutoff swept on |shrunk_logFC|.
# Source data: the benchmark cutoff_sweep.csv (shrunk_absLFC rows = lfsr<0.05 & |shrunk_logFC|>cut),
# computed on the current canonical_deg_results.csv.

suppressMessages({library(data.table); library(ggplot2); library(patchwork)})

ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
theme_file <- file.path(ROOT, "scripts/figures/publication_theme.R")
if (file.exists(theme_file)) try(source(theme_file), silent = TRUE)
base_theme <- if (exists("theme_pub")) theme_pub() else
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
base_theme <- base_theme + theme(plot.title = element_text(face = "plain", colour = "black"),
                                 axis.text = element_text(colour = "black"),
                                 axis.title = element_text(colour = "black"),
                                 legend.text = element_text(colour = "black"),
                                 legend.title = element_text(colour = "black"))

sweep <- fread(file.path(ROOT,
  "RNA-seq/results/audit_sensitivity/lfc_cutoff_benchmark/cutoff_sweep.csv"))
d <- sweep[scale == "shrunk_absLFC"][order(threshold)]

OUT <- file.path(ROOT,
  "figures/supplementary/figS_methods_validation/lfc_sensitivity")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# reference cutoffs to annotate
refs <- c(0.3, 0.5)

# ---- Panel 1: DEG set size vs cutoff (log y) ----
p_size <- ggplot(d, aes(threshold, nDEG)) +
  geom_vline(xintercept = refs, linetype = "dashed", colour = "#9E9E9E") +
  geom_line(colour = "black") + geom_point(colour = "black", size = 1.8) +
  geom_text(aes(label = scales::comma(nDEG)), vjust = -0.8, size = 2.6, colour = "black") +
  scale_y_log10(labels = scales::comma, expand = expansion(mult = c(0.05, 0.12))) +
  scale_x_continuous(breaks = d$threshold) +
  labs(title = "DEG set size vs ashr-shrunk |log2FC| cutoff (lfsr < 0.05)",
       x = NULL, y = "N DEGs (log scale)") +
  base_theme

# ---- Panel 2: known-MASLD recall vs cutoff ----
rl <- melt(d, id.vars = "threshold",
           measure.vars = c("recall_known", "recall_govaere", "recall_feng", "recall_steato"),
           variable.name = "panel", value.name = "recall")
lab <- c(recall_known = "Known union", recall_govaere = "Govaere",
         recall_feng = "Feng", recall_steato = "SteatoSITE")
rl[, panel := factor(lab[as.character(panel)], levels = lab)]
pal <- c("Known union" = "#000000", "Govaere" = "#1B7837",
         "Feng" = "#2166AC", "SteatoSITE" = "#B2182B")

p_recall <- ggplot(rl, aes(threshold, recall, colour = panel)) +
  geom_vline(xintercept = refs, linetype = "dashed", colour = "#9E9E9E") +
  geom_line(aes(linewidth = panel == "Known union")) +
  geom_point(size = 1.6) +
  scale_linewidth_manual(values = c(`TRUE` = 1.1, `FALSE` = 0.5), guide = "none") +
  scale_colour_manual(values = pal, name = NULL) +
  scale_x_continuous(breaks = d$threshold) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0.02, 0.05))) +
  labs(title = "Known-MASLD recall vs ashr-shrunk |log2FC| cutoff",
       x = "|shrunk log2FC| cutoff", y = "Recall of known genes (%)") +
  base_theme + theme(legend.position = "right")

fig <- p_size / p_recall + plot_layout(heights = c(1, 1))
outpdf <- file.path(OUT, "shrunk_cutoff_decision.pdf")
ggsave(outpdf, fig, width = 8, height = 7, device = cairo_pdf)
fwrite(d, file.path(OUT, "shrunk_cutoff_decision_data.csv"))

# caption to stdout (not in figure)
message("CAPTION: Decision panel for the canonical ashr-shrunk |log2FC| floor. ",
        "Significance gate fixed at lfsr<0.05 (canonical DEG definition); the |shrunk_logFC| ",
        "cutoff is swept on x. Top: pooled-analysis DEG set size. Bottom: recall of curated ",
        "MASLD gene panels (Govaere, Feng, SteatoSITE) and their union, among DEGs at each cutoff. ",
        "Dashed grey lines mark 0.3 and 0.5 for reference. Source: cutoff_sweep.csv (shrunk_absLFC), ",
        "computed on canonical_deg_results.csv.")
message("Wrote: ", outpdf)
print(d[, .(threshold, nDEG, recall_known, recall_govaere, recall_feng, recall_steato)])
