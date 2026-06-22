#!/usr/bin/env Rscript
# Sister to fig3_intro_cascade.pdf — colocalization tier expressed as #
# GWAS loci with >=1 COLOC-positive gene (PP4>=0.5), instead of total
# colocalized genes. Same canonical definitions and palette.
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

BASE <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
OUT  <- file.path(BASE, "figures/sketches/fig3_intro")
DAT  <- file.path(BASE, "scripts/figures/sketches_fig3_intro/data")
source(file.path(BASE, "scripts/figures/publication_theme.R"))

cascade <- fread(file.path(DAT, "method_cascade_loci.csv"))
cascade[, step := factor(step, levels = step)]
cascade[, group := fcase(
  grepl("^Lead loci|SuSiE converged|SuSiE-X|meSuSiE", step), "Finemapping",
  grepl("^Loci",                                       step), "Colocalization")]
cascade[, group := factor(group, levels = c("Finemapping","Colocalization"))]

step_fill <- c(
  "Lead loci screened"                                = "#90CAF9",
  "SuSiE converged loci (canonical)"                  = "#42A5F5",
  "SuSiE-X high-PIP physical loci (cross-ancestry)"   = "#1976D2",
  "meSuSiE shared CS physical loci (cross-ancestry)"  = "#0D47A1",
  "Loci — regular (ABF) COLOC tested"             = "#F48FB1",
  "Loci — SuSiE COLOC tested"                     = "#C2185B")

p <- ggplot(cascade, aes(x = count, y = step, fill = step)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = scales::label_comma()(count)),
            hjust = -0.18, size = 3.8, fontface = "bold") +
  scale_x_continuous(labels = scales::label_comma(),
                     expand = expansion(mult = c(0, 0.18))) +
  scale_y_discrete(limits = rev) +
  scale_fill_manual(values = step_fill, guide = "none") +
  facet_grid(group ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(x = "Count", y = NULL,
       subtitle = "Lead loci → finemapping → cross-ancestry refinement → loci where COLOC was tested") +
  theme_masld() +
  theme(plot.subtitle = element_text(size = 11, color = "gray30"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 11),
        strip.text.y.left = element_text(angle = 0, face = "bold", size = 11),
        strip.placement = "outside",
        strip.background = element_blank(),
        panel.spacing.y = unit(0.4, "cm"))

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(OUT, "option_10b_cascade_loci.pdf"),
       p, width = 11, height = 5, device = cairo_pdf)

# Also drop into the main fig3 panels dir
FIG3_PANELS <- file.path(BASE, "figures/main/fig2_genetics/panels")
ggsave(file.path(FIG3_PANELS, "intro_cascade_loci.pdf"),
       p, width = 11, height = 5, device = cairo_pdf)

cat("Wrote: ", file.path(OUT, "option_10b_cascade_loci.pdf"),
    "\n   and: ", file.path(FIG3_PANELS, "intro_cascade_loci.pdf"), "\n")
