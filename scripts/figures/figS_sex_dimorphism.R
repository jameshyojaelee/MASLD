#!/usr/bin/env Rscript
##############################################################################
# figS_sex_dimorphism.R  (total v3 overhaul, 2026-05-14)
#
# KEY MESSAGE (biological): MASLD disease responses are overwhelmingly
# concordant between sexes (>99% of expressed genes). The handful of robust
# sex-specific disease responses are AUTOSOMAL — not chrX/chrY-driven — and
# F_only carries an unusual share of non-coding (lncRNA / pseudogene)
# transcripts, suggesting a regulatory rather than coding axis. Companion
# figure figS_coloc_sex_direct connects these RNA-level signatures to GWAS /
# COLOC evidence and shows the sex × causal-genetics intersection is empty.
#
# Eight panels (4 rows × 2 cols):
#   a  Sample size F vs M across cohorts × healthy/disease          (DESIGN)
#   b  mashr shrinkage map: raw β_int vs shrunken β_int             (METHOD)
#   c  β_F vs β_M landscape colored by v3 assigned_class            (BIOLOGY)
#   d  Chromosomal category × class                                 (BIOLOGY)
#   e  Top 8 F_only genes (β_F dominant)                            (BIOLOGY)
#   f  Top 8 divergent genes (opposing β_F, β_M)              (BIOLOGY)
#   g  Gene-biotype composition by class (lncRNA-rich F_only)       (BIOLOGY)
#   h  Power-matched bootstrap stability by class                  (METHOD)
##############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(ggrepel)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

FIGDIR <- FIGS_SEX_DIR
OUT    <- file.path(FIGDIR, "figS_sex_dimorphism.pdf")

# ---------------------------------------------------------------------------
# Inputs (v3 canonical + v2 baseline frozen 2026-05-13)
# ---------------------------------------------------------------------------
V3_CSV   <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration",
  "sex_v3/sex_deg_classification_v3.csv")
V3_INPUT <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration",
  "sex_v3/intermediates/sex_v3_input.rds")
stopifnot(file.exists(V3_CSV), file.exists(V3_INPUT))
v3 <- fread(V3_CSV)
v3_input <- readRDS(V3_INPUT)
info <- as.data.table(v3_input$meta)

# ---------------------------------------------------------------------------
# Palette: v3 5-class (triple-gate adds an explicit "uncertain" class)
# ---------------------------------------------------------------------------
v3_colors <- c(
  F_only     = "#AD1457",   # rose magenta (female)
  M_only     = "#1A237E",   # navy (male — empty in v3)
  divergent  = "#6A1B9A",   # purple
  concordant = "#9E9E9E",   # gray
  uncertain  = "#CFD8DC"    # light blue-gray (post-triple-gate residual)
)
v3_levels <- c("F_only", "M_only", "divergent", "concordant", "uncertain")
# (mashr native column posterior_P_anti_correlated retained for provenance;
#  assigned_class_gated is the post-hoc triple-gate label.)

# Class label vector for compact axis use
v3_labels <- c(F_only     = "F_only",
               M_only     = "M_only",
               divergent  = "divergent",
               concordant = "concordant",
               uncertain  = "uncertain")

# Use post-hoc gated class label (sign-concordance gate + anti->divergent rename)
v3[, assigned_class := factor(assigned_class_gated, levels = v3_levels)]

# Cohort short-name mapping (5-cohort canonical)
COHORT_NAMES <- c(
  GSE126848 = "Suppli",
  GSE130970 = "Hoang",
  GSE135251 = "Govaere",
  GSE162694 = "Bril",
  GSE213621 = "Chen"
)

# ---------------------------------------------------------------------------
# Panel A: Sample size F vs M across cohorts × healthy/disease (DESIGN)
# Motivates the shrinkage: 99 F controls vs 50 M controls (2:1 imbalance)
# ---------------------------------------------------------------------------
message("Panel A: sample sizes M vs F × healthy/disease ...")

info[, dataset := factor(dataset, levels = names(COHORT_NAMES))]
info[, cohort_short := COHORT_NAMES[as.character(dataset)]]
info[, sex_group := factor(paste(inferred_sex, group_binary, sep = " "),
                            levels = c("F Control", "F Disease",
                                       "M Control", "M Disease"))]
n_dat <- info[, .N, by = .(cohort_short, sex_group)]
# Order cohorts by total sample size
order_n <- info[, .N, by = cohort_short][order(-N), cohort_short]
n_dat[, cohort_short := factor(cohort_short, levels = order_n)]

# Compute totals for headline annotation
tot_FC <- info[inferred_sex == "F" & group_binary == "Control", .N]
tot_FD <- info[inferred_sex == "F" & group_binary == "Disease", .N]
tot_MC <- info[inferred_sex == "M" & group_binary == "Control", .N]
tot_MD <- info[inferred_sex == "M" & group_binary == "Disease", .N]
ratio_ctrl <- round(tot_FC / tot_MC, 2)

sg_pal <- c(
  "F Control" = "#F8BBD0",   # soft pink
  "F Disease" = "#AD1457",   # rose magenta
  "M Control" = "#90CAF9",   # soft blue
  "M Disease" = "#1A237E"    # navy
)

panel_A <- ggplot(n_dat, aes(x = cohort_short, y = N, fill = sex_group)) +
  geom_col(width = 0.65, color = "white", linewidth = 0.2) +
  geom_text(data = info[, .N, by = cohort_short],
            aes(x = cohort_short, y = N, label = paste0("n=", N)),
            inherit.aes = FALSE,
            vjust = -0.4, size = 2.0, fontface = "bold", color = "gray25") +
  scale_fill_manual(values = sg_pal, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title    = "5-cohort canonical sample sizes",
    subtitle = sprintf("Controls: %d F vs %d M (F:M = %.2f); Disease: %d F vs %d M",
                       tot_FC, tot_MC, ratio_ctrl, tot_FD, tot_MD),
    x = NULL, y = "Samples"
  ) +
  theme_masld(base_size = 7) +
  theme(legend.position = "top",
        legend.key.size = unit(8, "pt"),
        legend.text     = element_text(size = 6))

save_fig(panel_A, file.path(FIGDIR, "figS_sex_a_sample_size.pdf"),
         width = fig_half_width, height = 3)

# ---------------------------------------------------------------------------
# Panel B: mashr shrinkage map — raw β_int vs posterior-mean β_int (METHOD)
# Visualizes what mashr does: weak interactions get pulled toward y=0;
# strong, replicable interactions retain their magnitude.
# ---------------------------------------------------------------------------
message("Panel B: mashr shrinkage map ...")

shrink <- copy(v3)
shrink[, raw_int     := beta_F - beta_M]
shrink[, shrunken_int := posterior_mean_F - posterior_mean_M]

# Downsample bulk classes (uncertain, concordant) for visual breath
set.seed(42)
unc_idx  <- which(shrink$assigned_class == "uncertain")
conc_idx <- which(shrink$assigned_class == "concordant")
oth_idx  <- which(!shrink$assigned_class %in% c("uncertain", "concordant"))
keep_idx <- c(sample(unc_idx,  min(length(unc_idx),  3000)),
              sample(conc_idx, min(length(conc_idx), 3000)),
              oth_idx)
shrink_p <- shrink[keep_idx]

lim_b <- ceiling(quantile(abs(shrink$raw_int), probs = 0.999,
                          na.rm = TRUE) * 2) / 2

# Median shrinkage by class
med_raw_conc  <- median(abs(shrink[assigned_class == "concordant"]$raw_int),
                         na.rm = TRUE)
med_shr_conc  <- median(abs(shrink[assigned_class == "concordant"]$shrunken_int),
                         na.rm = TRUE)
med_raw_fonly <- median(abs(shrink[assigned_class == "F_only"]$raw_int),
                         na.rm = TRUE)
med_shr_fonly <- median(abs(shrink[assigned_class == "F_only"]$shrunken_int),
                         na.rm = TRUE)

panel_B <- ggplot(shrink_p,
                  aes(x = raw_int, y = shrunken_int,
                      color = assigned_class)) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "gray70") +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "gray70") +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.4, linetype = "dashed",
              color = "gray55") +
  # Layer bottom-to-top: uncertain (lightest, hidden), concordant (gray),
  # then F_only/divergent/M_only on top (vibrant colors visible)
  rasterize_layer(
    geom_point(data = shrink_p[assigned_class == "uncertain"],
               alpha = 0.12, size = 0.25),
    dpi = 300
  ) +
  rasterize_layer(
    geom_point(data = shrink_p[assigned_class == "concordant"],
               alpha = 0.2, size = 0.35),
    dpi = 300
  ) +
  rasterize_layer(
    geom_point(data = shrink_p[!assigned_class %in% c("uncertain", "concordant")],
               alpha = 0.9, size = 0.9),
    dpi = 300
  ) +
  scale_color_manual(values = v3_colors, drop = TRUE,
                     name = NULL, labels = v3_labels) +
  coord_cartesian(xlim = c(-lim_b, lim_b), ylim = c(-lim_b, lim_b)) +
  labs(
    title    = "mashr shrinkage map",
    subtitle = sprintf("median |β_int|: concordant %.2f→%.3f; F_only %.2f→%.2f",
                       med_raw_conc, med_shr_conc,
                       med_raw_fonly, med_shr_fonly),
    x = expression("Raw" ~ beta[F] - beta[M] ~ "(dream)"),
    y = expression("Shrunken" ~ beta[F] - beta[M] ~ "(mashr posterior)")
  ) +
  guides(color = guide_legend(override.aes = list(size = 1.8, alpha = 1))) +
  theme_masld(base_size = 7) +
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = alpha("white", 0.7),
                                         color = NA),
        legend.key.size = unit(8, "pt"),
        legend.text     = element_text(size = 6))

save_fig(panel_B, file.path(FIGDIR, "figS_sex_b_shrinkage_map.pdf"),
         width = fig_half_width, height = 3)

# ---------------------------------------------------------------------------
# Panel C: β_F vs β_M LANDSCAPE (the biological map)
# ---------------------------------------------------------------------------
message("Panel A: β_F vs β_M landscape ...")

set.seed(42)
unc_idx  <- which(v3$assigned_class == "uncertain")
conc_idx <- which(v3$assigned_class == "concordant")
oth_idx  <- which(!v3$assigned_class %in% c("uncertain", "concordant"))
keep_idx <- c(sample(unc_idx,  min(length(unc_idx),  4000)),
              sample(conc_idx, min(length(conc_idx), 4000)),
              oth_idx)
scat_p   <- v3[keep_idx]

lim_v <- ceiling(quantile(c(abs(v3$beta_F), abs(v3$beta_M)),
                          probs = 0.998, na.rm = TRUE) * 2) / 2

# Highlight a handful of biological landmarks for direct labeling
landmark_F <- v3[assigned_class == "F_only"][
  order(-posterior_P_F_only)][1:3, gene_symbol]
landmark_A <- v3[assigned_class == "divergent"][
  order(-posterior_P_anti_correlated)][1:3, gene_symbol]
landmarks <- v3[gene_symbol %in% c(landmark_F, landmark_A)]

panel_C <- ggplot(scat_p,
                  aes(x = beta_M, y = beta_F, color = assigned_class)) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "gray70") +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "gray70") +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.4, linetype = "dashed",
              color = "gray55") +
  # Layer bottom-to-top: uncertain (faint), concordant (gray), then sex-specific (vibrant)
  rasterize_layer(
    geom_point(data = scat_p[assigned_class == "uncertain"],
               alpha = 0.10, size = 0.25),
    dpi = 300
  ) +
  rasterize_layer(
    geom_point(data = scat_p[assigned_class == "concordant"],
               alpha = 0.20, size = 0.35),
    dpi = 300
  ) +
  rasterize_layer(
    geom_point(data = scat_p[!assigned_class %in% c("uncertain", "concordant")],
               alpha = 0.9, size = 0.9),
    dpi = 300
  ) +
  geom_text_repel(data = landmarks,
                  aes(label = gene_symbol), size = 2.0,
                  segment.size = 0.2, segment.color = "gray40",
                  min.segment.length = 0, max.overlaps = Inf,
                  box.padding = 0.4, point.padding = 0.2,
                  show.legend = FALSE) +
  scale_color_manual(values = v3_colors, drop = TRUE,
                     name = NULL, labels = v3_labels) +
  coord_cartesian(xlim = c(-lim_v, lim_v), ylim = c(-lim_v, lim_v)) +
  labs(
    title    = "Per-gene disease effect by sex",
    subtitle = sprintf("27,638 genes | dashed = equal-effect (β_F = β_M); 99.4%% sit near the diagonal"),
    x = expression(beta[M] ~ "(disease effect in males)"),
    y = expression(beta[F] ~ "(disease effect in females)")
  ) +
  guides(color = guide_legend(override.aes = list(size = 1.8, alpha = 1))) +
  theme_masld(base_size = 7) +
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = alpha("white", 0.7),
                                         color = NA),
        legend.key.size = unit(8, "pt"),
        legend.text     = element_text(size = 6))

save_fig(panel_C, file.path(FIGDIR, "figS_sex_c_beta_landscape.pdf"),
         width = fig_half_width, height = 3)

# ---------------------------------------------------------------------------
# Panel D: Chromosomal distribution by class (autosomal vs chrX vs chrY)
# ---------------------------------------------------------------------------
message("Panel D: chromosomal distribution by class ...")

v3[, chr_bin := fcase(
  chr == "chrY", "chrY",
  chr == "chrX", "chrX",
  chr == "chrM", "chrM",
  default       = "autosomal"
)]
chr_bin_levels <- c("autosomal", "chrX", "chrY", "chrM")
v3[, chr_bin := factor(chr_bin, levels = chr_bin_levels)]

b_dat <- v3[!is.na(assigned_class), .N, by = .(assigned_class, chr_bin)]
# Fill missing class × chr cells with 0 (v5 vocab: F_only, M_only, divergent,
# concordant, uncertain)
b_dat <- merge(
  CJ(assigned_class = c("F_only", "M_only", "divergent", "concordant", "uncertain"),
     chr_bin        = chr_bin_levels,
     unique = TRUE),
  b_dat, by = c("assigned_class", "chr_bin"), all.x = TRUE
)
b_dat[is.na(N), N := 0]
b_dat[, assigned_class := factor(assigned_class,
        levels = c("F_only", "M_only", "divergent", "concordant", "uncertain"))]
b_dat[, chr_bin := factor(chr_bin, levels = rev(chr_bin_levels))]

# Compute per-class totals for axis labels
class_totals <- v3[!is.na(assigned_class), .N, by = assigned_class]
get_n <- function(cl) {
  v <- class_totals[assigned_class == cl, N]
  if (length(v) == 0) 0L else v
}
chr_xlabs <- c(
  F_only     = sprintf("F_only\n(n=%s)",     format(get_n("F_only"),     big.mark = ",")),
  M_only     = sprintf("M_only\n(n=%s)",     format(get_n("M_only"),     big.mark = ",")),
  divergent  = sprintf("divergent\n(n=%s)",  format(get_n("divergent"),  big.mark = ",")),
  concordant = sprintf("concordant\n(n=%s)", format(get_n("concordant"), big.mark = ",")),
  uncertain  = sprintf("uncertain\n(n=%s)",  format(get_n("uncertain"),  big.mark = ","))
)

panel_D <- ggplot(b_dat, aes(x = assigned_class, y = chr_bin,
                              fill = log10(N + 1))) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = format(N, big.mark = ","),
                color = log10(N + 1) > 2.5),
            size = 2.4, fontface = "bold") +
  scale_fill_gradient(low = "#F5F5F5", high = "#37474F",
                      name = expression(log[10](count + 1))) +
  scale_color_manual(values = c("TRUE" = "white", "FALSE" = "gray20"),
                     guide = "none") +
  scale_x_discrete(position = "top", labels = chr_xlabs) +
  labs(
    title    = "Sex-specific responses are autosomal",
    subtitle = "Counts per class × chromosome — F_only carries no chrX/chrY/chrM genes",
    x = NULL, y = NULL
  ) +
  theme_masld(base_size = 7) +
  theme(legend.position = "right",
        legend.key.size = unit(6, "pt"),
        legend.title    = element_text(size = 6),
        legend.text     = element_text(size = 6),
        axis.text.x     = element_text(face = "bold"),
        axis.text.y     = element_text(face = "bold"),
        axis.ticks      = element_blank(),
        panel.grid      = element_blank())

save_fig(panel_D, file.path(FIGDIR, "figS_sex_d_chr_distribution.pdf"),
         width = fig_half_width, height = 3)

# ---------------------------------------------------------------------------
# Panel E: Top 8 F_only genes (β_F dominant biology)
# ---------------------------------------------------------------------------
message("Panel E: top F_only genes ...")

top_F <- v3[assigned_class == "F_only"][
  order(-posterior_P_F_only)][1:10,
  .(gene_symbol, chr, post_p = posterior_P_F_only, beta_F, beta_M)]
top_F[, gene_symbol := ifelse(is.na(gene_symbol) | gene_symbol == "",
                               sub("^ENSG", "ENSG…", gene_symbol),
                               gene_symbol)]
top_F[, gene_lab := sprintf("%s  P=%.2f", gene_symbol, post_p)]
top_F <- unique(top_F, by = "gene_lab")
top_F[, gene_lab := factor(gene_lab, levels = rev(unique(gene_lab)))]

topFl <- melt(top_F, id.vars = c("gene_lab"),
              measure.vars = c("beta_F", "beta_M"),
              variable.name = "sex_eff", value.name = "beta")
topFl[, sex_eff := factor(sex_eff,
                          levels = c("beta_F", "beta_M"),
                          labels = c("β_F (female)", "β_M (male)"))]

panel_E <- ggplot() +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "gray60") +
  geom_segment(data = top_F,
               aes(y = gene_lab, yend = gene_lab,
                   x = beta_F, xend = beta_M),
               color = "#AD1457", linewidth = 0.5, alpha = 0.45) +
  geom_point(data = topFl,
             aes(y = gene_lab, x = beta, fill = sex_eff, shape = sex_eff),
             size = 2.0, color = "black", stroke = 0.25) +
  scale_fill_manual(values = c("β_F (female)" = "#AD1457",
                                "β_M (male)"   = "#1A237E"),
                    name = NULL) +
  scale_shape_manual(values = c("β_F (female)" = 21, "β_M (male)" = 24),
                     name = NULL) +
  labs(
    title    = "Top 10 F_only genes",
    x = expression(beta ~ "(disease vs control LFC)"),
    y = NULL
  ) +
  theme_masld(base_size = 7) +
  theme(axis.text.y     = element_text(size = 6),
        legend.position = "bottom",
        legend.key.size = unit(8, "pt"),
        legend.text     = element_text(size = 6))

save_fig(panel_E, file.path(FIGDIR, "figS_sex_e_top_fonly.pdf"),
         width = fig_half_width, height = 3)

# ---------------------------------------------------------------------------
# Panel F: Top 8 divergent genes (opposing-direction biology)
# ---------------------------------------------------------------------------
message("Panel F: top divergent genes ...")

top_A <- v3[assigned_class == "divergent"][
  order(-posterior_P_anti_correlated)][seq_len(min(10, .N)),
  .(gene_symbol, chr, post_p = posterior_P_anti_correlated, beta_F, beta_M)]
top_A <- top_A[!is.na(post_p)]
top_A[, gene_symbol := ifelse(is.na(gene_symbol) | gene_symbol == "",
                               sub("^ENSG", "ENSG…", gene_symbol),
                               gene_symbol)]
top_A[, gene_lab := sprintf("%s  P=%.2f", gene_symbol, post_p)]
top_A <- unique(top_A, by = "gene_lab")
top_A[, gene_lab := factor(gene_lab, levels = rev(unique(gene_lab)))]

topAl <- melt(top_A, id.vars = c("gene_lab"),
              measure.vars = c("beta_F", "beta_M"),
              variable.name = "sex_eff", value.name = "beta")
topAl[, sex_eff := factor(sex_eff,
                          levels = c("beta_F", "beta_M"),
                          labels = c("β_F (female)", "β_M (male)"))]

panel_F <- ggplot() +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "gray60") +
  geom_segment(data = top_A,
               aes(y = gene_lab, yend = gene_lab,
                   x = beta_F, xend = beta_M),
               color = "#6A1B9A", linewidth = 0.5, alpha = 0.45) +
  geom_point(data = topAl,
             aes(y = gene_lab, x = beta, fill = sex_eff, shape = sex_eff),
             size = 2.0, color = "black", stroke = 0.25) +
  scale_fill_manual(values = c("β_F (female)" = "#AD1457",
                                "β_M (male)"   = "#1A237E"),
                    name = NULL) +
  scale_shape_manual(values = c("β_F (female)" = 21, "β_M (male)" = 24),
                     name = NULL) +
  labs(
    title    = "Top 10 divergent genes",
    x = expression(beta ~ "(disease vs control LFC)"),
    y = NULL
  ) +
  theme_masld(base_size = 7) +
  theme(axis.text.y     = element_text(size = 6),
        legend.position = "bottom",
        legend.key.size = unit(8, "pt"),
        legend.text     = element_text(size = 6))

save_fig(panel_F, file.path(FIGDIR, "figS_sex_f_top_divergent.pdf"),
         width = fig_half_width, height = 3)

# ---------------------------------------------------------------------------
# Panel G: Gene-biotype composition by class
#   Biological insight: F_only and divergent classes carry an UNUSUAL
#   share of non-coding transcripts (lncRNA / pseudogene), distinguishing
#   them from the bulk-protein-coding concordant signal.
# ---------------------------------------------------------------------------
message("Panel G: gene_biotype composition by class ...")

# Collapse biotypes into 4 buckets for legibility
v3[, biotype_bin := fcase(
  gene_biotype == "protein_coding",                       "protein_coding",
  grepl("lncRNA|lincRNA|antisense", gene_biotype),         "lncRNA",
  grepl("pseudogene", gene_biotype),                       "pseudogene",
  default                                                  = "other"
)]
biotype_levels <- c("protein_coding", "lncRNA", "pseudogene", "other")
v3[, biotype_bin := factor(biotype_bin, levels = biotype_levels)]

e_dat <- v3[!is.na(assigned_class), .N, by = .(assigned_class, biotype_bin)]
e_dat[, total := sum(N), by = assigned_class]
e_dat[, frac  := N / total]
e_dat[, assigned_class := factor(assigned_class,
                                  levels = c("F_only", "M_only", "divergent",
                                             "concordant", "uncertain"))]

biotype_pal <- c(
  "protein_coding" = "#37474F",
  "lncRNA"         = "#AD1457",
  "pseudogene"     = "#FF9800",
  "other"          = "#BDBDBD"
)

panel_G <- ggplot(e_dat, aes(x = assigned_class, y = frac, fill = biotype_bin)) +
  geom_col(width = 0.65, color = "white", linewidth = 0.2) +
  geom_text(data = unique(e_dat[, .(assigned_class, total)]),
            aes(x = assigned_class, y = 1.04,
                label = sprintf("n=%s", format(total, big.mark = ","))),
            inherit.aes = FALSE,
            size = 2.1, fontface = "bold", color = "gray25") +
  scale_fill_manual(values = biotype_pal, name = NULL) +
  scale_y_continuous(labels = label_percent(),
                     limits = c(0, 1.10),
                     breaks = c(0, 0.25, 0.5, 0.75, 1.0),
                     expand = expansion(mult = c(0, 0))) +
  scale_x_discrete(labels = c(F_only = "F_only",
                              M_only = "M_only",
                              divergent = "divergent",
                              concordant = "concordant",
                              uncertain = "uncertain")) +
  labs(
    title    = "Non-coding RNAs over-represented in F_only",
    subtitle = "Gene-biotype composition by class",
    x = NULL, y = "Fraction of class"
  ) +
  theme_masld(base_size = 7) +
  theme(legend.position = "top",
        legend.key.size = unit(8, "pt"),
        legend.text     = element_text(size = 6))

save_fig(panel_G, file.path(FIGDIR, "figS_sex_g_biotype.pdf"),
         width = fig_half_width, height = 3)

# ---------------------------------------------------------------------------
# Panel H: METHOD HIGHLIGHT — Power-matched bootstrap stability per class
# Fraction of bootstrap reps in which each gene's TRIPLE-GATED class equals
# its full-data call. The pm_frac_* columns come from 04d aggregation of a
# 100-rep power-matched bootstrap (each rep subsamples the larger sex group
# to match the smaller before refitting dream + mashr + triple gate).
# Stacked bars per class show the distribution of stability bins.
# ---------------------------------------------------------------------------
message("Panel H: power-matched bootstrap stability per class ...")

## v5 100-rep bootstrap stability from bootstrap_stability_v5.csv.
## Per-class self-stability = stability_<class>_v5 for genes called as <class>.
boot_csv <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration",
  "sex_v3/bootstrap_stability_v5.csv")
have_pm <- file.exists(boot_csv)

if (have_pm) {
  boot <- fread(boot_csv,
                select = c("gene", "stability_Female_biased",
                           "stability_Male_biased", "stability_Divergent",
                           "stability_Divergent_one_sided",
                           "stability_Concordant", "stability_Sex_modifier",
                           "stability_Not_DEG"))
  v3 <- merge(v3, boot, by = "gene", all.x = TRUE)
  # Self-stability per v5 vocabulary
  v3[, pm_frac_self := fcase(
    assigned_class == "F_only",     as.numeric(stability_Female_biased),
    assigned_class == "M_only",     as.numeric(stability_Male_biased),
    assigned_class == "divergent",  as.numeric(pmax(stability_Divergent,
                                                    stability_Divergent_one_sided,
                                                    na.rm = TRUE)),
    assigned_class == "concordant", as.numeric(pmax(stability_Concordant,
                                                    stability_Sex_modifier,
                                                    na.rm = TRUE)),
    assigned_class == "uncertain",  as.numeric(stability_Not_DEG),
    default = NA_real_
  )]

  stab_levels <- c("0-25%", "25-50%", "50-80%", "80-100%")
  pm_dat <- v3[assigned_class %in% c("F_only", "M_only", "divergent",
                                      "concordant") &
                 !is.na(pm_frac_self)]
  pm_dat[, stab_bin := cut(pm_frac_self,
                            breaks = c(-Inf, 0.25, 0.50, 0.80, Inf),
                            labels = stab_levels,
                            right  = TRUE,
                            include.lowest = TRUE)]
  pm_dat[, assigned_class := factor(assigned_class,
                                     levels = c("F_only", "M_only", "divergent",
                                                "concordant"))]
  pm_dat[, stab_bin := factor(stab_bin, levels = stab_levels)]

  pm_counts <- pm_dat[, .N, by = .(assigned_class, stab_bin)]
  pm_counts[, total := sum(N), by = assigned_class]
  pm_counts[, frac  := N / total]

  stab_pal <- c(
    "0-25%"   = "#C62828",   # red (unstable)
    "25-50%"  = "#EF6C00",   # orange
    "50-80%"  = "#F9A825",   # amber
    "80-100%" = "#2E7D32"    # green (stable)
  )

  panel_H <- ggplot(pm_counts,
                    aes(x = assigned_class, y = frac, fill = stab_bin)) +
    geom_col(width = 0.65, color = "white", linewidth = 0.2) +
    geom_text(data = unique(pm_counts[, .(assigned_class, total)]),
              aes(x = assigned_class, y = 1.04,
                  label = sprintf("n=%s", format(total, big.mark = ","))),
              inherit.aes = FALSE,
              size = 2.1, fontface = "bold", color = "gray25") +
    scale_fill_manual(values = stab_pal,
                      name = "PM bootstrap stability\n(fraction reps same class)",
                      drop = FALSE) +
    scale_y_continuous(labels = label_percent(),
                       limits = c(0, 1.10),
                       breaks = c(0, 0.25, 0.5, 0.75, 1.0),
                       expand = expansion(mult = c(0, 0))) +
    labs(
      title    = "Power-matched bootstrap stability per class",
      subtitle = "Stacked share of genes by stability bin (100 PM reps)",
      x = NULL, y = "Fraction of class"
    ) +
    theme_masld(base_size = 7) +
    theme(legend.position = "top",
          legend.key.size = unit(7, "pt"),
          legend.text     = element_text(size = 6),
          legend.title    = element_text(size = 6))
} else {
  # Graceful placeholder if 04d aggregation hasn't been merged yet.
  panel_H <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = "pm_frac_* columns not yet present\n(run 04d_power_match_aggregate.R)",
             size = 2.6, color = "gray35") +
    theme_void()
}

save_fig(panel_H, file.path(FIGDIR, "figS_sex_h_pm_stability.pdf"),
         width = fig_half_width, height = 3)

# ---------------------------------------------------------------------------
# Composite (2 rows × 3 cols)
# ---------------------------------------------------------------------------
message("Composing figS_sex_dimorphism.pdf ...")

composite <- (panel_A | panel_B) /
             (panel_C | panel_D) /
             (panel_E | panel_F) /
             (panel_G | panel_H) +
  plot_annotation(
    title    = "Sex-dimorphic disease responses in MASLD: design, biology, statistical robustness",
    subtitle = "Bayesian sex classification (mashr on dream M2; cohort RE + multiply-imputed age + SVA + composition) | triple gate: posterior_P > 0.8, lfsr < 0.05, sign-concordance | 100-rep power-matched bootstrap | 794 samples, 27,638 genes",
    tag_levels = "a",
    theme = theme(
      plot.title    = element_text(size = 9,   face = "bold"),
      plot.subtitle = element_text(size = 6.8, color = "gray35"),
      plot.tag      = element_text(size = 8,   face = "bold")
    )
  ) &
  theme_masld(base_size = 7)

save_fig_tall(composite, OUT,
              width  = fig_full_width,
              height = 12)

message("=== Done ===")
message("Output: ", FIGDIR)
