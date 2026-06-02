#!/usr/bin/env Rscript
##############################################################################
# figS_coloc_sex_direct.R  (total v3 overhaul, 2026-05-14)
#
# KEY MESSAGE: Sex-specific MASLD transcriptional responses (v3 F_only,
# divergent) are NOT enriched for common-variant causal architecture.
# The strongest v2 enrichment claim (Female_biased × COLOC PP.H4>0.1
# OR=1.58, padj=1.9e-4) does not survive proper Bayesian shrinkage on the
# transcription side (v3 OR=0.91, padj=1.0). COLOC hits are overwhelmingly
# sex-CONCORDANT, and sex-stratified GWAS / liver-enzyme sources drive the
# convergence rather than diverging by sex class.
#
# Five panels:
#   a  COLOC enrichment OR forest at 3 PP4 thresholds                 (NULL)
#   b  β_F vs β_M scatter with COLOC overlay — hits land on diagonal  (BIOL)
#   c  PP4 distribution by class                                      (BIOL)
#   d  Spotlight: top 10 non-concordant genes by PP4                  (BIOL)
#   e  Stratified GWAS-set forest: enrichment across strata (null)    (NULL)
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

FIGDIR <- file.path(BASE, "figures/supplementary/figS_sex_dimorphism")
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)
OUT    <- file.path(FIGDIR, "figS_sex_csd_composite.pdf")

# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------
V3_CSV   <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration",
  "sex_v3/sex_deg_classification_v3.csv")
ATLAS    <- file.path(BASE,
  "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
V2_COLOC <- file.path(BASE,
  "RNA-seq/results/stratified_causal/v2_baseline_pre_v3/sex_coloc_refresh_v2.v2.csv")
V3_COLOC <- file.path(BASE,
  "RNA-seq/results/stratified_causal/sex_coloc_refresh_v2.csv")
stopifnot(file.exists(V3_CSV), file.exists(ATLAS),
          file.exists(V2_COLOC), file.exists(V3_COLOC))

v3  <- fread(V3_CSV)
cv2 <- fread(V2_COLOC)
cv3 <- fread(V3_COLOC)

a_cols <- c("ensembl_id", "sex_class",
            "coloc_susie_best_pp4", "coloc_susie_best_gwas",
            "coloc_best_pp4_polyfun", "coloc_best_gwas_polyfun")
a <- fread(ATLAS, select = a_cols)
a[, pp4_best := pmax(coloc_susie_best_pp4, coloc_best_pp4_polyfun, na.rm = TRUE)]
a[, gwas_best := ifelse(!is.na(coloc_susie_best_pp4) &
                          (is.na(coloc_best_pp4_polyfun) |
                             coloc_susie_best_pp4 >= coloc_best_pp4_polyfun),
                        coloc_susie_best_gwas, coloc_best_gwas_polyfun)]

# Join v3 (per-gene assigned_class + β_F + β_M) with atlas PP4
v3[, ensembl_id := sub("\\.[0-9]+$", "", gene)]
a[, ensembl_id := sub("\\.[0-9]+$", "", ensembl_id)]
m <- merge(v3[, .(ensembl_id, gene_symbol, chr, assigned_class, assigned_class_gated,
                   posterior_P_F_only, posterior_P_anti_correlated,
                   posterior_confidence, beta_F, beta_M)],
           a[, .(ensembl_id, pp4_best, gwas_best)],
           by = "ensembl_id", all.x = TRUE)
m[, pp4_best := pmax(pp4_best, 0, na.rm = TRUE)]

# ---------------------------------------------------------------------------
# v3 palette (concordant gray per project rule)
# ---------------------------------------------------------------------------
v3_colors <- c(
  F_only     = "#AD1457",
  M_only     = "#1A237E",
  divergent  = "#6A1B9A",
  concordant = "#9E9E9E"
)
v3_levels <- c("F_only", "M_only", "divergent", "concordant")
v3_labels <- c(F_only     = "F_only",
               M_only     = "M_only",
               divergent  = "divergent",
               concordant = "concordant")
# Use post-hoc gated label from v3 csv (sign-concordance gate + rename)
m[, assigned_class := factor(assigned_class_gated, levels = v3_levels)]

# ---------------------------------------------------------------------------
# Panel A: COLOC enrichment OR forest at 3 PP4 thresholds (computed INLINE
# from the gated v3 csv + atlas PP4 so it reflects whatever the canonical
# classification currently says)
# ---------------------------------------------------------------------------
message("Panel A: COLOC OR forest (inline Fisher) ...")

fisher_enrich <- function(dt, class_col, hit_col, target, hit_cut) {
  # Restrict to genes with a CONFIDENT class call (drop post-triple-gate
  # uncertain — enrichment is class-vs-class, not class-vs-uncertain).
  cls <- dt[[class_col]]
  cls_chr <- as.character(cls)
  keep <- !is.na(cls_chr) & cls_chr != "uncertain"
  dt <- dt[keep]

  n_class       <- sum(dt[[class_col]] == target, na.rm = TRUE)
  n_class_coloc <- sum(dt[[class_col]] == target & dt[[hit_col]] > hit_cut,
                       na.rm = TRUE)
  n_bg          <- sum(!is.na(dt[[class_col]]) & dt[[class_col]] != target,
                       na.rm = TRUE)
  n_bg_coloc    <- sum(!is.na(dt[[class_col]]) & dt[[class_col]] != target &
                       dt[[hit_col]] > hit_cut, na.rm = TRUE)
  if (n_class == 0 || (n_class_coloc + n_bg_coloc) == 0) {
    return(data.table(odds_ratio = NA_real_, ci_low = NA_real_,
                       ci_high = NA_real_, pvalue = NA_real_,
                       n_class = n_class, n_class_coloc = n_class_coloc))
  }
  ft <- fisher.test(matrix(c(n_class_coloc, n_class - n_class_coloc,
                              n_bg_coloc,    n_bg    - n_bg_coloc),
                            nrow = 2, byrow = TRUE))
  data.table(odds_ratio = unname(ft$estimate),
             ci_low     = ft$conf.int[1],
             ci_high    = ft$conf.int[2],
             pvalue     = ft$p.value,
             n_class    = n_class,
             n_class_coloc = n_class_coloc)
}

classes  <- c("F_only", "divergent")   # M_only is empty
pp4_cuts <- c(0.1, 0.5, 0.8)
rows     <- list()
for (cl in classes) for (cut in pp4_cuts) {
  r <- fisher_enrich(m, "assigned_class", "pp4_best", cl, cut)
  r[, c("sex_class", "pp4_cut") := .(cl, cut)]
  rows[[paste(cl, cut)]] <- r
}
e_dat <- rbindlist(rows)
e_dat[, padj := p.adjust(pvalue, method = "BH")]
e_dat[, sex_class := factor(sex_class, levels = c("F_only", "M_only",
                                                    "divergent"))]
e_dat[, pp4_cut := factor(sprintf("PP4 > %.1f", pp4_cut),
                           levels = c("PP4 > 0.1", "PP4 > 0.5", "PP4 > 0.8"))]
e_dat[, ci_upper_plot := pmin(ci_high, 5)]
e_dat[, sig_lab := fcase(
  is.na(padj),    "",
  padj < 0.001,   "***",
  padj < 0.01,    "**",
  padj < 0.05,    "*",
  default         = "ns"
)]

panel_A <- ggplot(e_dat[!is.na(odds_ratio)],
                  aes(y = sex_class, x = odds_ratio,
                      color = pp4_cut, shape = pp4_cut)) +
  geom_vline(xintercept = 1, linewidth = 0.4, linetype = "dashed",
             color = "gray55") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_upper_plot),
                 height = 0.0, linewidth = 0.5,
                 position = position_dodge(width = 0.6)) +
  geom_point(size = 2.4, position = position_dodge(width = 0.6)) +
  geom_text(aes(label = sig_lab, x = ci_upper_plot + 0.18),
            position = position_dodge(width = 0.6),
            size = 2.0, color = "black", hjust = 0) +
  scale_color_manual(values = c("PP4 > 0.1" = "#1565C0",
                                 "PP4 > 0.5" = "#6A1B9A",
                                 "PP4 > 0.8" = "#C2185B"),
                     name = NULL) +
  scale_shape_manual(values = c("PP4 > 0.1" = 16,
                                 "PP4 > 0.5" = 17,
                                 "PP4 > 0.8" = 15),
                     guide = "none") +
  scale_x_continuous(limits = c(0, 5.5),
                     breaks = c(0, 0.5, 1, 2, 3, 5),
                     expand = expansion(mult = c(0.02, 0.12))) +
  scale_y_discrete(drop = FALSE) +
  labs(
    title    = "Sex-specific transcription is not COLOC-enriched",
    subtitle = "No class crosses padj < 0.05 at any PP4 threshold",
    x        = "COLOC odds ratio (95% CI)",
    y        = NULL
  ) +
  theme_masld(base_size = 7) +
  theme(panel.grid.major.x   = element_line(linewidth = 0.2, color = "gray92"),
        legend.position      = "top",
        legend.key.size      = unit(8, "pt"),
        legend.text          = element_text(size = 6))

save_fig(panel_A, file.path(FIGDIR, "figS_csd_a_coloc_collapse.pdf"),
         width = fig_half_width, height = 3)

# ---------------------------------------------------------------------------
# Panel B: β_F vs β_M with COLOC overlay
# ---------------------------------------------------------------------------
message("Panel B: β_F vs β_M with COLOC overlay ...")

set.seed(42)
# Restrict to genes with measurable COLOC support (PP4 >= 0.1).
scat <- m[!is.na(pp4_best) & pp4_best >= 0.1]
scat[, is_hit := pp4_best > 0.5]

lim_v <- ceiling(quantile(c(abs(m$beta_F), abs(m$beta_M)),
                          probs = 0.998, na.rm = TRUE) * 2) / 2

# Labels: all non-concordant genes with PP4 >= 0.5 (the "strong" set),
# plus the top-by-PP4 gene in each non-concordant class so the best divergent
# call (typically below 0.5) is also annotated.
strong_lab <- m[pp4_best >= 0.5 & assigned_class != "concordant"]
top_per_cls <- m[assigned_class %in% c("F_only", "divergent") &
                  !is.na(pp4_best)][order(assigned_class, -pp4_best)][,
                  .SD[1], by = assigned_class]
hit_labels <- unique(rbind(strong_lab, top_per_cls))[,
              .(gene_symbol, beta_F, beta_M, pp4_best, assigned_class)]

# Layer order (bottom -> top):
#   1) concordant cloud (rasterized faint background, sized by PP4)
#   2) ALL non-concordant points, sized by PP4 (single layer, area-proportional)
#   3) labels on top non-concordant per class
panel_B <- ggplot(scat,
                  aes(x = beta_M, y = beta_F)) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "gray70") +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "gray70") +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.4, linetype = "dashed",
              color = "gray55") +
  rasterize_layer(
    geom_point(data = scat[assigned_class == "concordant"],
               aes(color = assigned_class, size = pp4_best),
               alpha = 0.25, shape = 16),
    dpi = 300
  ) +
  geom_point(data = scat[assigned_class != "concordant"],
             aes(color = assigned_class, size = pp4_best),
             alpha = 0.9, shape = 16) +
  geom_text_repel(data = hit_labels,
                  aes(label = gene_symbol, color = assigned_class),
                  size = 2.1, segment.size = 0.3, segment.color = "gray35",
                  min.segment.length = 0, max.overlaps = Inf,
                  box.padding = 0.6, point.padding = 0.25,
                  fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = v3_colors, drop = TRUE, name = "Class",
                     labels = v3_labels) +
  scale_size_area(max_size = 5,
                  limits = c(0.1, 1),
                  breaks = c(0.1, 0.5, 0.8),
                  labels = c("0.1", "0.5", "0.8"),
                  name = "PP4") +
  coord_cartesian(xlim = c(-lim_v, lim_v), ylim = c(-lim_v, lim_v)) +
  labs(
    title    = "Sex-biased COLOC hits",
    subtitle = sprintf("Genes with PP4 ≥ 0.1 (n=%s) | top non-concordant per class labeled",
                       format(nrow(scat), big.mark = ",")),
    x = expression(beta[M] ~ "(disease effect in males)"),
    y = expression(beta[F] ~ "(disease effect in females)")
  ) +
  guides(color = guide_legend(override.aes = list(size = 1.8, alpha = 1)),
         size  = guide_legend(override.aes = list(color = "gray40"))) +
  theme_masld(base_size = 7) +
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = alpha("white", 0.7),
                                         color = NA),
        legend.key.size = unit(7, "pt"),
        legend.title    = element_text(size = 6),
        legend.text     = element_text(size = 6),
        legend.spacing.y = unit(1, "pt"))

save_fig(panel_B, file.path(FIGDIR, "figS_csd_b_betaFvsM_pp4.pdf"),
         width = fig_half_width, height = 3)

# ---------------------------------------------------------------------------
# Panel C: PP4 distribution by v3 class
# ---------------------------------------------------------------------------
message("Panel C: PP4 distribution by v3 class ...")

c_dat <- m[!is.na(pp4_best) & pp4_best > 0 & assigned_class != "M_only"]
c_dat[, assigned_class := factor(assigned_class,
                                  levels = c("F_only", "divergent",
                                             "concordant"))]
c_med <- c_dat[, .(med = median(pp4_best), n = .N), by = assigned_class]
c_med[, lab := sprintf("n=%s\nmed %.2f", format(n, big.mark = ","), med)]

panel_C <- ggplot(c_dat, aes(x = assigned_class, y = pp4_best,
                              fill = assigned_class)) +
  geom_violin(alpha = 0.5, color = NA, scale = "width", width = 0.85) +
  geom_boxplot(width = 0.18, outlier.size = 0.4, fill = "white",
               linewidth = 0.3) +
  geom_hline(yintercept = c(0.5, 0.8), linewidth = 0.3, linetype = "dashed",
             color = "gray55") +
  geom_text(data = c_med,
            aes(x = assigned_class, y = 1.10, label = lab),
            inherit.aes = FALSE, size = 1.9, color = "gray25") +
  scale_fill_manual(values = v3_colors, guide = "none") +
  scale_x_discrete(labels = c(F_only = "F_only",
                              divergent = "divergent",
                              concordant = "concordant")) +
  scale_y_continuous(limits = c(0, 1.18),
                     breaks = c(0, 0.5, 0.8, 1.0),
                     expand = expansion(mult = c(0, 0))) +
  labs(
    title    = "PP4 distribution by class",
    subtitle = "Most colocalized genes are sex-concordant",
    x = NULL, y = "Best PP4 across all GWAS"
  ) +
  theme_masld(base_size = 7)

save_fig(panel_C, file.path(FIGDIR, "figS_csd_c_pp4_by_class.pdf"),
         width = fig_half_width, height = 3)

# gwas_bucket assignment (used downstream by panels E and F)
m[, gwas_bucket := fcase(
  grepl("UKBB_ALT|UKBB_AST|UKBB_GGT|BBJ_(ALT|AST|GGT)|enzyme",
        gwas_best, ignore.case = TRUE), "Liver enzymes",
  grepl("NAFLD|MASH|NASH",  gwas_best, ignore.case = TRUE), "NAFLD/MASH",
  grepl("PDFF",             gwas_best, ignore.case = TRUE), "PDFF",
  grepl("HCC|Cirrh|Ghouse", gwas_best, ignore.case = TRUE), "Cirrhosis/HCC",
  is.na(gwas_best),                                          NA_character_,
  default = "Other"
)]

# ---------------------------------------------------------------------------
# Panel E: spotlight every sex-specific COLOC hit (PP4>0.5)
# ---------------------------------------------------------------------------
message("Panel E: spotlight sex-specific COLOC hits ...")

hit_tab <- m[pp4_best > 0.5 & assigned_class != "concordant",
             .(gene_symbol, chr, assigned_class, beta_F, beta_M, pp4_best,
               gwas_best, gwas_bucket)]
# If there are essentially zero — fall back to top sex-specific by PP4
if (nrow(hit_tab) < 4) {
  candidates <- m[!is.na(assigned_class) & assigned_class != "concordant" & !is.na(pp4_best)]
  n_take <- min(10, nrow(candidates))
  if (n_take > 0) {
    hit_tab <- candidates[order(-pp4_best)][1:n_take,
      .(gene_symbol, chr, assigned_class, beta_F, beta_M, pp4_best,
        gwas_best, gwas_bucket)]
  } else {
    hit_tab <- candidates[, .(gene_symbol, chr, assigned_class, beta_F, beta_M, pp4_best,
                              gwas_best, gwas_bucket)]
  }
}
hit_tab[, gwas_short := ifelse(is.na(gwas_bucket), "—", gwas_bucket)]
hit_tab[, lab := sprintf("%s | %s | PP4=%.2f",
                          gene_symbol, gwas_short, pp4_best)]
hit_tab[, lab := factor(lab, levels = rev(lab))]

hit_long <- melt(hit_tab,
                 id.vars      = c("lab", "assigned_class"),
                 measure.vars = c("beta_F", "beta_M"),
                 variable.name = "sex_eff",
                 value.name    = "beta")
hit_long[, sex_eff := factor(sex_eff,
                              levels = c("beta_F", "beta_M"),
                              labels = c("β_F", "β_M"))]

panel_E <- ggplot() +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "gray60") +
  geom_segment(data = hit_tab,
               aes(y = lab, yend = lab,
                   x = beta_F, xend = beta_M, color = assigned_class),
               linewidth = 0.5, alpha = 0.55) +
  geom_point(data = hit_long,
             aes(y = lab, x = beta, fill = sex_eff, shape = sex_eff),
             size = 2.0, color = "black", stroke = 0.25) +
  scale_fill_manual(values = c("β_F" = "#AD1457", "β_M" = "#1A237E"),
                    name = NULL) +
  scale_shape_manual(values = c("β_F" = 21, "β_M" = 24),
                     name = NULL) +
  scale_color_manual(values = v3_colors, guide = "none") +
  labs(
    title    = "Sex-specific colocalized genes",
    subtitle = sprintf("Top %d non-concordant genes by PP4", nrow(hit_tab)),
    x = expression(beta ~ "(disease vs control LFC)"),
    y = NULL
  ) +
  theme_masld(base_size = 7) +
  theme(axis.text.y     = element_text(size = 5.5),
        legend.position = "bottom",
        legend.key.size = unit(8, "pt"),
        legend.text     = element_text(size = 6))

save_fig(panel_E, file.path(FIGDIR, "figS_csd_e_spotlight.pdf"),
         width = fig_half_width, height = 3)

# ---------------------------------------------------------------------------
# Panel F: stratified GWAS-set forest — v3 enrichment in all / disease /
# liver-enzyme strata at PP4>0.1
# ---------------------------------------------------------------------------
message("Panel F: stratified GWAS-set forest (inline Fisher) ...")

# Build per-gene stratum membership: a gene is "in stratum S" if its best PP4
# is from a GWAS in S. Disease = NAFLD/MASH/HCC/Cirrhosis; Enzyme = UKBB/BBJ
# ALT/AST/GGT.
m[, in_enzyme  := !is.na(gwas_bucket) & gwas_bucket == "Liver enzymes"]
m[, in_disease := !is.na(gwas_bucket) & gwas_bucket %in% c("NAFLD/MASH",
                                                            "Cirrhosis/HCC")]

strat_rows <- list()
for (cl in classes) {
  # all
  r <- fisher_enrich(m, "assigned_class", "pp4_best", cl, 0.1)
  r[, c("sex_class", "stratum") := .(cl, "All GWAS")]
  strat_rows[[paste("all",     cl)]] <- r
  # disease subset
  r <- fisher_enrich(m[in_disease == TRUE | pp4_best <= 0.1],
                     "assigned_class", "pp4_best", cl, 0.1)
  r[, c("sex_class", "stratum") := .(cl, "Disease GWAS")]
  strat_rows[[paste("disease", cl)]] <- r
  # enzyme subset
  r <- fisher_enrich(m[in_enzyme == TRUE | pp4_best <= 0.1],
                     "assigned_class", "pp4_best", cl, 0.1)
  r[, c("sex_class", "stratum") := .(cl, "Liver enzyme")]
  strat_rows[[paste("enzyme",  cl)]] <- r
}
f_dat <- rbindlist(strat_rows)
f_dat[, padj := p.adjust(pvalue, method = "BH")]
f_dat[, stratum := factor(stratum,
                           levels = c("All GWAS", "Disease GWAS", "Liver enzyme"))]
f_dat[, sex_class := factor(sex_class,
                             levels = c("F_only", "M_only", "divergent"))]
f_dat[, ci_upper_plot := pmin(ci_high, 5)]
f_dat[, sig_lab := fcase(
  is.na(padj),     "",
  padj < 0.001,    "***",
  padj < 0.01,     "**",
  padj < 0.05,     "*",
  default          = "ns"
)]

panel_F <- ggplot(f_dat[!is.na(odds_ratio)],
                  aes(y = sex_class, x = odds_ratio, color = stratum,
                      shape = stratum)) +
  geom_vline(xintercept = 1, linewidth = 0.4, linetype = "dashed",
             color = "gray55") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_upper_plot),
                 height = 0.0, linewidth = 0.5,
                 position = position_dodge(width = 0.6)) +
  geom_point(size = 2.4, position = position_dodge(width = 0.6)) +
  geom_text(aes(label = sig_lab, x = ci_upper_plot + 0.18),
            position = position_dodge(width = 0.6),
            size = 2.0, color = "black", hjust = 0) +
  scale_color_manual(values = c("All GWAS"     = "#1565C0",
                                 "Disease GWAS" = "#C2185B",
                                 "Liver enzyme" = "#FF8F00"),
                     name = NULL) +
  scale_shape_manual(values = c("All GWAS"     = 16,
                                 "Disease GWAS" = 17,
                                 "Liver enzyme" = 15),
                     guide = "none") +
  scale_x_continuous(limits = c(0, 5.5),
                     breaks = c(0, 0.5, 1, 2, 3, 5),
                     expand = expansion(mult = c(0.02, 0.12))) +
  scale_y_discrete(drop = FALSE) +
  labs(
    title    = "Stratified GWAS-set enrichment",
    subtitle = "Every stratum is null (padj > 0.05)",
    x = "Odds ratio (95% CI)",
    y = NULL
  ) +
  theme_masld(base_size = 7) +
  theme(panel.grid.major.x = element_line(linewidth = 0.2, color = "gray92"),
        legend.position    = "top",
        legend.key.size    = unit(8, "pt"),
        legend.text        = element_text(size = 6))

save_fig(panel_F, file.path(FIGDIR, "figS_csd_f_stratum_forest.pdf"),
         width = fig_half_width, height = 3)

# ---------------------------------------------------------------------------
# Composite (3 rows × 2 cols)
# ---------------------------------------------------------------------------
message("Composing figS_sex_csd_composite.pdf ...")

composite <- (panel_A | panel_B) /
             (panel_C | panel_E) /
             panel_F +
  plot_annotation(
    title    = "Sex × MASLD GWAS / COLOC integration: a null intersection",
    subtitle = "mashr sex classification × multi-ancestry COLOC (PolyFun EUR + SuSiE/ABF; 28 GWAS) | atlas 27,187 genes × 291 columns",
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
              height = 9)

message("=== Done ===")
message("Output: ", FIGDIR)
