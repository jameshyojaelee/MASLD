#!/usr/bin/env Rscript
# KEY MESSAGE: Currin 2025 138-donor bulk-liver caQTL independently validates
# motif-disrupting variants — direction-concordant in 80-83% of HNF4A/THRB/RORA
# variants, and 47 variants reach 4-way concordance (motif + eQTL + DA + caQTL).

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

OUT_PDF <- file.path(FIGS05_DIR, "figS05_scatac_caqtl_concordance.pdf")

GA_DIR <- file.path(BASE, "GWAS/finemapping/results/gwas_atac")
var_dt <- fread(file.path(GA_DIR, "currin_caqtl_concordance.csv"))
tf_dt  <- fread(file.path(GA_DIR, "currin_caqtl_per_tf.csv"))

# ---- Color palette (Liang style) -----------------------------------------
MAG  <- "#C9265E"
BLUE <- "#1565C0"
GOLD <- "#FFB300"
GRAY <- "#9E9E9E"

FOURWAY <- c("THRB", "HNF4A", "RORA", "MLXIPL", "MAX")

# ---- Panel a: variant-level motif x caQTL effect scatter -----------------
# Keep variants tested in caQTL with non-missing motif + caQTL effects
v <- var_dt[currin_caqtl_tested == TRUE & !is.na(alleleDiff) & !is.na(currin_caqtl_beta)]
v[, concord_class := fcase(
  concord_4of4 == TRUE,        "4-way concordant",
  concord_motif_caqtl == TRUE, "Motif x caQTL concordant",
  default = "Discordant"
)]
v[, concord_class := factor(concord_class,
  levels = c("Discordant", "Motif x caQTL concordant", "4-way concordant"))]
setorder(v, concord_class)  # plot 4-way on top

# Label the most extreme 4-way concordant variants by TF (one per top TF)
top_tf <- intersect(FOURWAY, v[concord_4of4 == TRUE]$tf_name)
label_dt <- v[concord_4of4 == TRUE & tf_name %in% top_tf,
              .SD[which.max(abs(alleleDiff) * abs(currin_caqtl_beta))],
              by = tf_name]

p_a <- ggplot(v, aes(x = alleleDiff, y = currin_caqtl_beta)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "gray70",
             linewidth = 0.25) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "gray70",
             linewidth = 0.25) +
  geom_point(aes(colour = concord_class), shape = 16,
             size = 0.7, alpha = 0.6) +
  geom_point(data = v[concord_4of4 == TRUE],
             aes(colour = concord_class), shape = 16, size = 1.4) +
  geom_text_repel(data = label_dt,
                  aes(label = tf_name),
                  size = 2.3, fontface = "bold", colour = "black",
                  segment.size = 0.2, segment.colour = "gray50",
                  min.segment.length = 0, box.padding = 0.3,
                  max.overlaps = Inf) +
  scale_colour_manual(values = c(
    "Discordant"              = GRAY,
    "Motif x caQTL concordant" = BLUE,
    "4-way concordant"         = MAG)) +
  scale_x_continuous(name = expression(bold("Motif disruption (alleleDiff)")),
                     limits = c(-1, 1) * max(abs(v$alleleDiff)),
                     expand = expansion(mult = 0.04)) +
  scale_y_continuous(name = expression(bold(paste("caQTL ", beta))),
                     expand = expansion(mult = 0.05)) +
  labs(tag = "a", colour = NULL) +
  theme_masld(base_size = 7) +
  theme_pub() +
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(0.25, "cm"),
        plot.tag = element_text(size = 11, face = "bold"),
        plot.tag.position = c(0.02, 0.97))

# ---- Panel b: per-TF concordance lollipop --------------------------------
tf <- tf_dt[n_caqtl_tested >= 4 & !is.na(pct_motif_caqtl_agree)]
# Always retain the 4-way headline set (even if MAX has n=1)
extra <- tf_dt[tf_name %in% FOURWAY & !(tf_name %in% tf$tf_name)]
tf <- rbind(tf, extra, fill = TRUE)
tf[, is_4way := tf_name %in% FOURWAY]
tf[, has_dz  := any_disease_regulon == TRUE]
# Rank: prefer disease-regulon + 4-way; sort by % concordance descending
setorder(tf, -is_4way, -has_dz, -pct_motif_caqtl_agree, -n_caqtl_tested)
tf_top <- head(tf, 18)
tf_top[, tf_name := factor(tf_name, levels = rev(tf_top$tf_name))]
tf_top[, col := fcase(is_4way == TRUE, MAG,
                     has_dz   == TRUE, BLUE,
                     default          = GRAY)]

p_b <- ggplot(tf_top, aes(y = tf_name)) +
  geom_vline(xintercept = 50, linetype = "dashed",
             colour = "gray60", linewidth = 0.3) +
  geom_segment(aes(x = 0, xend = pct_motif_caqtl_agree, yend = tf_name,
                   colour = col), linewidth = 0.45) +
  geom_point(data = tf_top[is_4way == TRUE],
             aes(x = pct_motif_caqtl_agree), colour = GOLD,
             size = 3.2, shape = 16) +
  geom_point(aes(x = pct_motif_caqtl_agree, size = n_caqtl_tested,
                 colour = col, fill = col), shape = 21, stroke = 0.5) +
  scale_colour_identity() +
  scale_fill_identity() +
  scale_size_continuous(range = c(1.2, 3.6),
                        breaks = c(5, 15, 40),
                        name = expression(bold(italic(n)*" tested"))) +
  scale_x_continuous(limits = c(0, 105),
                     breaks = c(0, 25, 50, 75, 100),
                     name = "caQTL concordance (%)",
                     expand = expansion(mult = c(0.01, 0.02))) +
  labs(tag = "b", y = NULL) +
  theme_masld(base_size = 7) +
  theme_pub() +
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(0.25, "cm"),
        axis.text.y = element_text(face = "bold", colour = "black"),
        axis.title.x = element_text(face = "bold"),
        plot.tag = element_text(size = 11, face = "bold"),
        plot.tag.position = c(0.02, 0.97))

# ---- Panel c: 4-way evidence breakdown ----------------------------------
# Per-TF n_4of4 count (TFs with at least one 4-way concordant variant)
fourway_tf <- tf_dt[n_4of4 >= 1]
setorder(fourway_tf, -n_4of4, -pct_motif_caqtl_agree)
fourway_tf <- head(fourway_tf, 15)
fourway_tf[, tf_name := factor(tf_name, levels = rev(fourway_tf$tf_name))]
fourway_tf[, is_4way := tf_name %in% FOURWAY]
fourway_tf[, has_dz  := any_disease_regulon == TRUE]
fourway_tf[, fill_col := fcase(is_4way == TRUE, MAG,
                               has_dz   == TRUE, BLUE,
                               default          = GRAY)]

p_c <- ggplot(fourway_tf, aes(y = tf_name, x = n_4of4)) +
  geom_col(aes(fill = fill_col), width = 0.65) +
  scale_fill_identity() +
  scale_x_continuous(name = expression(bold("4-way concordant variants  ("*italic(n)*")")),
                     expand = expansion(mult = c(0, 0.08)),
                     breaks = scales::breaks_pretty(n = 4)) +
  labs(tag = "c", y = NULL) +
  theme_masld(base_size = 7) +
  theme_pub() +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(face = "bold", colour = "black"),
        axis.title.x = element_text(face = "bold"),
        plot.tag = element_text(size = 11, face = "bold"),
        plot.tag.position = c(0.02, 0.97))

# ---- Compose -------------------------------------------------------------
combo <- p_a + p_b + p_c + plot_layout(widths = c(2.4, 2.6, 2.0))

ggsave(OUT_PDF, combo, width = 7.2, height = 3.0, device = cairo_pdf)
cat(sprintf("[saved] %s\n", OUT_PDF))
cat(sprintf("  variants tested in caQTL: %d (4-way: %d; motif x caQTL: %d)\n",
            nrow(v),
            sum(v$concord_4of4 == TRUE, na.rm = TRUE),
            sum(v$concord_motif_caqtl == TRUE, na.rm = TRUE)))
cat(sprintf("  TFs in panel b (top 18): %s\n",
            paste(rev(levels(tf_top$tf_name)), collapse = ", ")))
cat(sprintf("  TFs with >=1 4-way variant: %d total; top in panel c: %s\n",
            nrow(tf_dt[n_4of4 >= 1]),
            paste(rev(levels(fourway_tf$tf_name)), collapse = ", ")))
