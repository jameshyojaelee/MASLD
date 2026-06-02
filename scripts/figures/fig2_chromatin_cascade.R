#!/usr/bin/env Rscript
# KEY MESSAGE: Chromatin remodeling is late-stage and hepatocyte-dominant.
# Dots per (cell type, F-stage transition) with direct CT labels.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT_PDF  <- file.path(FIGS05_DIR, "figS05_scatac_chromatin_cascade.pdf")
DATA_DIR <- file.path(FIGS05_DIR, "data")
dir.create(dirname(OUT_PDF), recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)

STAGE_DA <- file.path(BASE, "Analysis/ATAC/Human_Multiome/results/stage_da")

cell_types  <- c("Hep", "Mac", "Fib", "Endo", "Chol")
transitions <- c("F0_vs_F2F3", "F0_vs_F4", "F2F3_vs_F4")
trans_lab   <- c(F0_vs_F2F3 = "F0→F2/F3",
                 F0_vs_F4   = "F0→F4",
                 F2F3_vs_F4 = "F2/F3→F4")
ct_label    <- c(Hep = "Hepatocytes", Mac = "Macrophages",
                 Fib = "Fibroblasts", Endo = "Endothelial",
                 Chol = "Cholangiocytes")
ct_pal      <- c(Hepatocytes    = "#C9265E",
                 Macrophages    = "#9C27B0",
                 Fibroblasts    = "#FF9800",
                 Endothelial    = "#1565C0",
                 Cholangiocytes = "#4CAF50")

count_da <- function(ct, tr) {
  f <- file.path(STAGE_DA, sprintf("da_%s_%s.csv", tr, ct))
  if (!file.exists(f)) return(NULL)
  dt <- fread(f, select = c("padj"))
  data.table(cell_type = ct_label[ct],
             transition = tr,
             n_sig = sum(dt$padj < 0.05, na.rm = TRUE))
}

dat <- rbindlist(lapply(transitions, function(tr)
  rbindlist(lapply(cell_types, function(ct) count_da(ct, tr)))))
dat[, transition := factor(trans_lab[transition], levels = unname(trans_lab))]
dat[, cell_type  := factor(cell_type, levels = unname(ct_label))]
fwrite(dat, file.path(DATA_DIR, "fig2_chromatin_cascade_main.csv"))

# Right-edge CT labels: place at the rightmost transition
label_dat <- dat[transition == levels(dat$transition)[3]]

panel <- ggplot(dat, aes(x = transition, y = n_sig + 1,
                         colour = cell_type, group = cell_type)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 2, shape = 16) +
  geom_text_repel(data = label_dat, aes(label = cell_type),
                  nudge_x = 0.25, hjust = 0, direction = "y",
                  segment.size = 0.2, segment.colour = "grey60",
                  size = 2.3, fontface = "bold", show.legend = FALSE,
                  min.segment.length = 0) +
  scale_colour_manual(values = ct_pal, guide = "none") +
  scale_y_log10(breaks = c(1, 10, 100, 1000),
                labels = c("0", "10", "100", "1000")) +
  scale_x_discrete(expand = expansion(add = c(0.4, 1.2))) +
  labs(x = NULL, y = "Significant DA peaks") +
  theme_masld(base_size = 7) +
  theme_pub() +
  theme(axis.title.y = element_text(face = "bold"),
        axis.text.x  = element_text(face = "bold"),
        panel.grid   = element_blank())

ggsave(OUT_PDF, panel, width = 3.6, height = 2.4, device = cairo_pdf)
cat(sprintf("[saved] %s\n", OUT_PDF))
