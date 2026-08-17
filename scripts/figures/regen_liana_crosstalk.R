#!/usr/bin/env Rscript
# Regenerate just fig2j (LIANA crosstalk panel)

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

fig2_data <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/fig2_data")
panel_dir <- file.path(FIG4_SC_DIR, "panels", "supplementary")
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)

liana_links_file <- file.path(fig2_data, "liana_ligand_tf_links.csv")
stopifnot(file.exists(liana_links_file))

links <- fread(liana_links_file)
links <- links[!is.na(link_score)]

links_dedup <- links[, .(
  link_score        = link_score[1],
  ligand_score_diff = ligand_score_diff[1],
  n_target_TFs      = uniqueN(target_TF),
  target_TFs         = paste(unique(target_TF), collapse = "/")
), by = .(sender_cell_type, ligand, receptor)]

top_links <- links_dedup[order(-link_score)][1:min(15, .N)]
top_links[, sender_short := sub(" cells$| derived cells$", "", sender_cell_type)]
top_links[, pair_label := paste0(ligand, " \u2192 ", receptor, " (", sender_short, ")")]
top_links[, pair_label := factor(pair_label, levels = rev(unique(top_links$pair_label)))]

sender_colors <- c("Cholangiocytes" = "#CE93D8", "Fibroblasts" = "#FFB74D",
                    "Macrophages" = "#E57373", "Endothelial" = "#80CBC4")

p_j <- ggplot(top_links,
              aes(x = abs(ligand_score_diff), y = pair_label,
                  fill = sender_short)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = sender_colors, name = "Sender") +
  theme_masld() +
  theme(axis.text.y  = element_text(size = 5, face = "italic"),
        axis.text.x  = element_text(size = 4.5),
        legend.position = "bottom",
        legend.key.size = unit(0.15, "cm"),
        legend.text     = element_text(size = 4.5),
        legend.margin = margin(0, 0, 0, 0),
        plot.margin = margin(5, 5, 5, 5)) +
  guides(fill = guide_legend(nrow = 2)) +
  labs(x = "|Ligand effect| (MASLD vs Ctrl)",
       y = NULL,
       title = "Hepatocyte-targeted L\u2013R\ninteractions (LIANA)")

out_file <- file.path(panel_dir, "figs3_liana_crosstalk.pdf")
ggsave(out_file, p_j, width = 4, height = 3.5, device = cairo_pdf)
message("Saved: ", out_file)
