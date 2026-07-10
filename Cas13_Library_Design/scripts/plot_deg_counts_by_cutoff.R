#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FIGLIB <- file.path(ROOT, "Cas13_Library_Design/figures")
DATDIR <- file.path(ROOT, "Cas13_Library_Design/data")
dir.create(FIGLIB, showWarnings = FALSE, recursive = TRUE)

RAW <- "raw:  padj<0.05 & |log2FC| > x"
ASHR <- "ashr: lfsr<0.05 & |shrunk| > x"
TREAT <- "treat: FDR<0.05 with lfc = x"

base_counts <- fread(file.path(DATDIR, "patient_concordance_by_cutoff_to3.csv"))[
  grepl("^Disease \\(", group) & (grepl("^raw:", condition) | grepl("^ashr:", condition)),
  .(condition, cutoff, n_degs)
]
base_counts[grepl("^raw:", condition), condition := RAW]
base_counts[grepl("^ashr:", condition), condition := ASHR]
treat_counts <- fread(file.path(DATDIR, "treat_deg_counts.csv"))[
  ,
  .(condition = TREAT, cutoff = treat_lfc, n_degs = total)
]

ndeg <- rbind(base_counts, treat_counts, fill = TRUE)
ndeg[, condition := factor(condition, levels = c(RAW, ASHR, TREAT))]
pal <- setNames(c("#B2182B", "#3B4CC0", "#1B9E77"), c(RAW, ASHR, TREAT))
bt <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = "top",
    legend.box = "vertical",
    plot.title = element_text(face = "bold", size = 11)
  )

p <- ggplot(ndeg, aes(cutoff, n_degs, color = condition)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.0) +
  scale_color_manual(values = pal) +
  scale_x_continuous(breaks = seq(0, 3, 0.5), limits = c(0, 3)) +
  scale_y_continuous(labels = comma) +
  labs(
    x = "Effect-size cutoff",
    y = "Total DEGs (up + down)",
    title = "Total DEG count by cutoff"
  ) +
  bt

out <- file.path(FIGLIB, "deg_counts_by_cutoff.pdf")
ggsave(out, p, width = 6.8, height = 4.2, useDingbats = FALSE)
cat("Wrote", out, "\n")
