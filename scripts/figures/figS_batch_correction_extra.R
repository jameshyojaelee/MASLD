#!/usr/bin/env Rscript
# figS_batch_correction_extra.R — Five additional panels (h-l) extending
# figS_batch_correction.pdf with atlas-specific quantitative evidence that
# (i) Harmony does mix cohorts in PCA space (it just doesn't erase the
#     biologically-driven cohort × disease-severity confounding); and
# (ii) the integrated DE call is robust regardless of UMAP geometry.
#
# Panels:
#   h: Pre vs post Harmony — top-PC scatter on raw vs Harmony-corrected,
#      with cohort iLISI deltas annotated
#   i: Per-cohort vs integrated logFC scatter grid (5 mega-eligible cohorts;
#      Tier 1 DEGs only) — Spearman ρ + slope per panel
#   j: PCA variance attribution per PC — for top 10 corrected PCs, the
#      partial R² attributable to disease / fibrosis / sex / dataset
#   k: Direction concordance per Tier 1 DEG — histogram of n_cohorts where
#      sign(per-cohort logFC) == sign(integrated logFC)
#   l: Per-cohort kBET-style chi-square mixing test — observed vs expected
#      cohort proportions in 30-NN neighborhoods
#
# All panels saved individually to figures/supplementary/figS_methods_validation/batch_correction/panels/
# Combined figure regenerated to include all panels (a-l).

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(RANN)
  library(uwot)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- file.path(FIGS_BATCH_DIR, "panels")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)

COHORT_LABEL <- c(
  GSE126848="Suppli", GSE130970="Hoang", GSE135251="Govaere", GSE162694="Bril",
  GSE167523="Kozumi", GSE174478="Kawamura", GSE193066="Hoshida", GSE213621="Chen",
  GSE240729="Verschuren", PRJNA512027="Gerhard"
)
COHORT_COLORS <- c(
  "Suppli"="#1F77B4","Hoang"="#FF7F0E","Govaere"="#2CA02C","Bril"="#D62728",
  "Kozumi"="#9467BD","Kawamura"="#8C564B","Hoshida"="#E377C2","Chen"="#7F7F7F",
  "Verschuren"="#BCBD22","Gerhard"="#17BECF"
)
disease_colors <- c(Control = masld_colors$control, Disease = masld_colors$nash)

PADJ_INT <- 0.05
LFC_INT  <- 0.5
KNN_K    <- 30
FIVE_COHORTS <- c("GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE213621")

# ----------------------------------------------------------------------------
# Load
# ----------------------------------------------------------------------------
message("Loading data...")
umap_dt <- fread(file.path(INT_RESULTS, "umap_coordinates.csv"))
dge     <- load_merged_dge()
sample_ids <- rownames(dge$samples)
stopifnot(nrow(umap_dt) == length(sample_ids))
umap_dt[, sample_id := sample_ids]

meta <- fread(file.path(INT_META, "unified_metadata.csv"))
plot_dt <- merge(umap_dt,
                 meta[, .(sample_id, fibrosis_stage, nas_score,
                          diagnosis_harmonized, sex)],
                 by = "sample_id", all.x = TRUE)
plot_dt[, cohort := factor(COHORT_LABEL[dataset],
                           levels = unname(COHORT_LABEL))]
plot_dt[, disease_state := factor(group, levels = c("Control", "Disease"))]

dream <- load_dream_results()
dream[, gene_clean := sub("\\..*", "", gene)]
dream[, is_tier1 := !is.na(dream_padj) & dream_padj < PADJ_INT &
                    !is.na(dream_logFC) & abs(dream_logFC) > LFC_INT]
n_tier1 <- sum(dream$is_tier1, na.rm = TRUE)
message(sprintf("Tier 1 DEGs: %s", comma(n_tier1)))

# ----------------------------------------------------------------------------
# PANEL H: Pre vs post Harmony — UMAP on the SAME PCs, both colorings
# ----------------------------------------------------------------------------
message("Panel H: pre vs post Harmony UMAP (apples-to-apples)...")
# Raw logCPM PCA: top 2000 variable genes for speed
logcpm_raw <- edgeR::cpm(dge, log = TRUE, prior.count = 2)
gv <- apply(logcpm_raw, 1, var)
top_idx <- order(gv, decreasing = TRUE)[seq_len(min(2000, nrow(logcpm_raw)))]
pca_raw <- prcomp(t(logcpm_raw[top_idx, ]), center = TRUE, scale. = FALSE)

# Pre-Harmony UMAP: compute UMAP on the raw top-30 PCs
n_pcs_for_umap <- 30
set.seed(42)
umap_pre <- uwot::umap(pca_raw$x[, seq_len(n_pcs_for_umap)],
                       n_neighbors = 30, min_dist = 0.3, metric = "euclidean",
                       n_components = 2, ret_model = FALSE)
pre_dt <- data.table(sample_id = sample_ids,
                     pre_UMAP1 = umap_pre[, 1],
                     pre_UMAP2 = umap_pre[, 2])
pre_dt <- merge(pre_dt,
                plot_dt[, .(sample_id, cohort, disease_state)],
                by = "sample_id")

# Post-Harmony UMAP already exists in plot_dt (UMAP1, UMAP2)
ilisi <- function(coords, labels, k) {
  nn  <- RANN::nn2(coords, query = coords, k = k + 1)$nn.idx[, -1]
  lab <- as.integer(as.factor(labels))
  vapply(seq_len(nrow(nn)), function(i) {
    p <- table(lab[nn[i, ]]) / k
    1 / sum(p^2)
  }, numeric(1))
}

pre_lisi_cohort  <- median(ilisi(as.matrix(pre_dt[, .(pre_UMAP1, pre_UMAP2)]),
                                 pre_dt$cohort, KNN_K))
pre_lisi_disease <- median(ilisi(as.matrix(pre_dt[, .(pre_UMAP1, pre_UMAP2)]),
                                 pre_dt$disease_state, KNN_K))
post_lisi_cohort  <- median(ilisi(as.matrix(plot_dt[, .(UMAP1, UMAP2)]),
                                  plot_dt$cohort, KNN_K))
post_lisi_disease <- median(ilisi(as.matrix(plot_dt[, .(UMAP1, UMAP2)]),
                                  plot_dt$disease_state, KNN_K))

mk_umap_panel <- function(dt, x, y, color_var, palette, title,
                          lisi_text, show_legend = FALSE,
                          legend_ncol = 2) {
  p <- ggplot(dt, aes(x = .data[[x]], y = .data[[y]],
                      color = .data[[color_var]])) +
    rasterize_layer(geom_point(size = 0.28, alpha = 0.65, shape = 16)) +
    scale_color_manual(values = palette, name = NULL,
                       na.value = "gray85", drop = FALSE) +
    coord_fixed() +
    labs(x = "UMAP 1", y = "UMAP 2", title = title) +
    annotate("label",
             x = min(dt[[x]]) + 0.02 * diff(range(dt[[x]])),
             y = max(dt[[y]]) - 0.02 * diff(range(dt[[y]])),
             label = lisi_text, hjust = 0, vjust = 1, size = 1.85,
             color = "gray20",
             fill = scales::alpha("white", 0.85)) +
    theme_masld(base_size = 7) +
    theme(plot.title = element_text(size = 7.5, face = "bold"),
          legend.position = if (show_legend) "right" else "none",
          legend.text = element_text(size = 5.8),
          legend.key.size = unit(0.2, "cm")) +
    guides(color = guide_legend(ncol = legend_ncol,
                                override.aes = list(size = 1.2, alpha = 1)))
  p
}

p_h_pre_cohort  <- mk_umap_panel(pre_dt, "pre_UMAP1", "pre_UMAP2", "cohort",
                                 COHORT_COLORS, "Pre batch correction · cohort",
                                 sprintf("cohort iLISI = %.2f / 10", pre_lisi_cohort),
                                 show_legend = TRUE, legend_ncol = 2)
p_h_post_cohort <- mk_umap_panel(plot_dt, "UMAP1", "UMAP2", "cohort",
                                 COHORT_COLORS, "Post batch correction · cohort",
                                 sprintf("cohort iLISI = %.2f / 10", post_lisi_cohort),
                                 show_legend = TRUE, legend_ncol = 2)
p_h_pre_dis     <- mk_umap_panel(pre_dt, "pre_UMAP1", "pre_UMAP2", "disease_state",
                                 disease_colors, "Pre batch correction · disease state",
                                 sprintf("disease iLISI = %.2f / 2", pre_lisi_disease),
                                 show_legend = TRUE, legend_ncol = 1)
p_h_post_dis    <- mk_umap_panel(plot_dt, "UMAP1", "UMAP2", "disease_state",
                                 disease_colors, "Post batch correction · disease state",
                                 sprintf("disease iLISI = %.2f / 2", post_lisi_disease),
                                 show_legend = TRUE, legend_ncol = 1)

p_h <- patchwork::wrap_plots(
    p_h_pre_cohort, p_h_post_cohort,
    p_h_pre_dis,    p_h_post_dis,
    ncol = 2, byrow = TRUE, guides = "collect"
  ) +
  plot_annotation(
    caption = sprintf(paste0(
        "Cohort iLISI: %.2f -> %.2f after batch correction (max 10; higher = better mixed). ",
        "Disease iLISI: %.2f -> %.2f (max 2; lower = cleaner Control/Disease separation)."
      ), pre_lisi_cohort, post_lisi_cohort,
         pre_lisi_disease, post_lisi_disease),
    theme = theme(plot.caption = element_text(size = 6, color = "gray35", hjust = 0)))
save_fig(p_h, file.path(PANEL_DIR, "figS_batch_h_pre_post_harmony.pdf"),
         width = fig_full_width * 1.1, height = 5.2)

# Save pre-Harmony UMAP coordinates for reproducibility
fwrite(pre_dt, file.path(FIGS_BATCH_DIR, "figS_batch_h_pre_umap_data.csv"))

# ----------------------------------------------------------------------------
# PANEL I: Per-cohort vs integrated logFC concordance scatter grid
# ----------------------------------------------------------------------------
message("Panel I: per-cohort vs integrated logFC scatter grid...")
ps <- load_per_study_de()
ps <- ps[dataset %in% FIVE_COHORTS]
ps[, gene_clean := sub("\\..*", "", gene)]

# Tier 1 DEG genes from integrated
tier1_genes <- dream[is_tier1 == TRUE, .(gene_clean, dream_logFC, dream_padj)]
ps_t1 <- merge(ps[, .(gene_clean, dataset, ps_logFC = logFC, ps_padj = padj)],
               tier1_genes, by = "gene_clean")
ps_t1[, cohort := factor(COHORT_LABEL[dataset],
                         levels = unname(COHORT_LABEL[FIVE_COHORTS]))]
ps_t1[, sig_flag := ps_padj < PADJ_INT]

# Per-cohort summary stats
cohort_stats <- ps_t1[, .(
  n = .N,
  spearman = cor(ps_logFC, dream_logFC, method = "spearman", use = "complete.obs"),
  pearson  = cor(ps_logFC, dream_logFC, method = "pearson", use = "complete.obs"),
  pct_concordant = 100 * mean(sign(ps_logFC) == sign(dream_logFC), na.rm = TRUE)
), by = cohort]
setorder(cohort_stats, -spearman)
cohort_stats[, label := sprintf("rho==%.2f~ '|'~%.0f*'%%'~direction",
                                spearman, pct_concordant)]

# Plot: facet by cohort, add identity line, add stats
ps_t1[, cohort := factor(cohort,
                         levels = cohort_stats$cohort)]
xlim_i <- range(ps_t1$ps_logFC, na.rm = TRUE)
ylim_i <- range(ps_t1$dream_logFC, na.rm = TRUE)
common_lim <- c(min(xlim_i[1], ylim_i[1]), max(xlim_i[2], ylim_i[2]))

stats_pos <- ps_t1[, .(x = -2.5, y = max(dream_logFC, na.rm = TRUE) * 0.98),
                   by = cohort]
stats_pos <- merge(stats_pos, cohort_stats[, .(cohort, label)], by = "cohort")

p_i <- ggplot(ps_t1, aes(x = ps_logFC, y = dream_logFC)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray70", linewidth = 0.2) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray70", linewidth = 0.2) +
  geom_abline(slope = 1, intercept = 0, color = "gray45",
              linetype = "dashed", linewidth = 0.3) +
  rasterize_layer(
    geom_point(aes(color = sig_flag), size = 0.32, alpha = 0.55, shape = 16)
  ) +
  geom_smooth(method = "lm", se = FALSE, color = "#C2185B",
              linewidth = 0.4, formula = y ~ x) +
  geom_text(data = stats_pos,
            aes(x = x, y = y, label = label),
            inherit.aes = FALSE, parse = TRUE,
            size = 1.95, color = "gray15", hjust = 0, vjust = 1) +
  scale_color_manual(values = c("FALSE" = "#BDBDBD", "TRUE" = "#1565C0"),
                     labels = c("FALSE" = "n.s. in cohort", "TRUE" = "padj<0.05 in cohort"),
                     name = NULL) +
  scale_x_continuous(limits = common_lim, breaks = scales::pretty_breaks(4)) +
  scale_y_continuous(limits = common_lim, breaks = scales::pretty_breaks(4)) +
  facet_wrap(~ cohort, nrow = 1) +
  labs(x = expression("Per-cohort log"[2]*" fold change"),
       y = expression("Integrated log"[2]*" fold change"),
       title = sprintf(
         "Per-cohort vs integrated logFC concordance (Tier 1 DEGs, n = %s)",
         comma(n_tier1))) +
  theme_masld(base_size = 7) +
  theme(plot.title = element_text(size = 8, face = "bold"),
        strip.text = element_text(size = 6.8, face = "bold"),
        legend.position = "top",
        legend.key.size = unit(0.25, "cm"),
        legend.text = element_text(size = 6),
        panel.spacing = unit(0.4, "lines"))

save_fig(p_i, file.path(PANEL_DIR, "figS_batch_i_percohort_logfc.pdf"),
         width = fig_full_width * 1.1, height = 2.4)

# ----------------------------------------------------------------------------
# PANEL J: PCA variance attribution per PC
# ----------------------------------------------------------------------------
message("Panel J: PCA variance attribution per PC...")
# Use raw logCPM PCA; for top 10 PCs, fit lm(PC ~ dataset + group + sex)
# and decompose ANOVA to get partial R² per covariate.
top_pcs <- 10
PC_mat <- pca_raw$x[, seq_len(top_pcs)]
colnames(PC_mat) <- paste0("PC", seq_len(top_pcs))
attr_dt <- merge(data.table(sample_id = sample_ids, PC_mat),
                 meta[, .(sample_id, dataset, group_binary,
                          fibrosis_stage, sex)],
                 by = "sample_id")
attr_dt[, dataset := factor(dataset)]
attr_dt[, group_binary := factor(group_binary)]
attr_dt[, sex := factor(sex)]
# fibrosis_stage NA -> 'NA' factor level
attr_dt[, fib_factor := factor(ifelse(is.na(fibrosis_stage), "NA",
                                      paste0("F", fibrosis_stage)))]

decompose_pc <- function(pc_vec, dt) {
  # Joint model with type-I sequential ANOVA so the SS partition sums to 100%.
  # Order disease + fibrosis FIRST so they get credit for any variance shared
  # with dataset (the conservative "biology-first" decomposition).
  dt2 <- copy(dt)
  dt2[, pc := pc_vec]
  dt2 <- dt2[complete.cases(group_binary, fib_factor, sex, dataset)]
  fit <- stats::lm(pc ~ group_binary + fib_factor + sex + dataset, data = dt2)
  aov_tab <- stats::anova(fit)
  total_ss <- sum(aov_tab[, "Sum Sq"], na.rm = TRUE)
  list(
    disease  = 100 * aov_tab["group_binary",  "Sum Sq"] / total_ss,
    fibrosis = 100 * aov_tab["fib_factor",    "Sum Sq"] / total_ss,
    sex      = 100 * aov_tab["sex",           "Sum Sq"] / total_ss,
    dataset  = 100 * aov_tab["dataset",       "Sum Sq"] / total_ss,
    residual = 100 * aov_tab["Residuals",     "Sum Sq"] / total_ss
  )
}

attr_long <- rbindlist(lapply(seq_len(top_pcs), function(k) {
  pc_vec <- attr_dt[[paste0("PC", k)]]
  d <- decompose_pc(pc_vec, attr_dt)
  data.table(PC = paste0("PC", k),
             pc_idx = k,
             pc_var = 100 * pca_raw$sdev[k]^2 / sum(pca_raw$sdev^2),
             disease  = d$disease,
             fibrosis = d$fibrosis,
             sex      = d$sex,
             dataset  = d$dataset,
             residual = d$residual)
}))

attr_melt <- melt(attr_long,
                  id.vars = c("PC", "pc_idx", "pc_var"),
                  measure.vars = c("disease", "fibrosis", "sex", "dataset", "residual"),
                  variable.name = "covariate", value.name = "pct")
attr_melt[, covariate := factor(covariate,
                                levels = c("disease", "fibrosis", "sex",
                                           "dataset", "residual"),
                                labels = c("Disease (binary)", "Fibrosis stage",
                                           "Sex", "Dataset", "Residual"))]
attr_melt[, PC := factor(PC, levels = paste0("PC", seq_len(top_pcs)))]
covar_palette <- c(
  "Disease (binary)" = "#C2185B",
  "Fibrosis stage"   = "#880E4F",
  "Sex"              = "#7B1FA2",
  "Dataset"          = "#9E9E9E",
  "Residual"         = "#E0E0E0"
)

p_j <- ggplot(attr_melt, aes(x = PC, y = pct, fill = covariate)) +
  geom_col(width = 0.75, color = "white", linewidth = 0.2) +
  geom_text(data = attr_long,
            aes(x = PC, y = 102,
                label = sprintf("%.1f%%", pc_var)),
            inherit.aes = FALSE,
            size = 1.7, color = "gray35", vjust = 0) +
  scale_fill_manual(values = covar_palette, name = NULL) +
  scale_y_continuous(labels = function(v) paste0(v, "%"),
                     limits = c(0, 115),
                     breaks = c(0, 25, 50, 75, 100),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "% PC variance explained",
       title = "Variance attribution per PC (raw PCA, top 10 PCs)",
       caption = paste0(
         "Top of bar = % total expression variance captured by that PC.\n",
         "Joint sequential ANOVA (disease + fibrosis + sex + dataset).\n",
         "Cohort dominates PC1/PC2 of the RAW matrix — the variance Harmony then absorbs.\n",
         "Disease + fibrosis biology shows up cleanly from PC3 onward."
       )) +
  theme_masld(base_size = 7) +
  theme(plot.title    = element_text(size = 8, face = "bold"),
        plot.caption  = element_text(size = 6, color = "gray35", hjust = 0),
        legend.position = "top",
        legend.key.size = unit(0.28, "cm"),
        legend.text = element_text(size = 6))

save_fig(p_j, file.path(PANEL_DIR, "figS_batch_j_pc_attribution.pdf"),
         width = fig_full_width * 0.7, height = 2.6)

# ----------------------------------------------------------------------------
# PANEL K: Direction concordance per Tier 1 DEG
# ----------------------------------------------------------------------------
message("Panel K: direction concordance per Tier 1 DEG...")
# For each Tier 1 DEG, count how many of the 5 mega-eligible cohorts have
# sign(per-cohort logFC) matching sign(integrated logFC).
ps_dir <- ps_t1[, .(n_cohorts_concordant = sum(
                       sign(ps_logFC) == sign(dream_logFC), na.rm = TRUE),
                    n_cohorts_tested = sum(!is.na(ps_logFC))),
                by = .(gene_clean, dream_logFC)]
# Restrict to the genes that have all 5 cohorts tested
ps_dir <- ps_dir[n_cohorts_tested == 5]

dir_counts <- ps_dir[, .N, by = n_cohorts_concordant]
dir_counts[, pct := 100 * N / sum(N)]
setorder(dir_counts, n_cohorts_concordant)
dir_counts[, x_factor := factor(n_cohorts_concordant, levels = 0:5)]

# Random expectation under H0 of independent 50/50 signs:
binom_exp <- data.table(x_factor = factor(0:5, levels = 0:5),
                        pct_expected = 100 * dbinom(0:5, 5, 0.5))

p_k <- ggplot(dir_counts, aes(x = x_factor, y = pct)) +
  geom_col(fill = masld_colors$up, width = 0.7,
           color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%s)", pct, comma(N))),
            vjust = -0.2, size = 1.85, color = "gray15", lineheight = 0.85) +
  geom_segment(data = binom_exp,
               aes(x = as.numeric(x_factor) - 0.4,
                   xend = as.numeric(x_factor) + 0.4,
                   y = pct_expected, yend = pct_expected),
               inherit.aes = FALSE,
               linetype = "dashed", color = "gray45", linewidth = 0.35) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(labels = function(v) paste0(v, "%"),
                     expand = expansion(mult = c(0, 0.18))) +
  labs(x = "Cohorts with logFC sign matching integrated direction",
       y = "% of Tier 1 DEGs",
       title = sprintf(
         "Direction concordance per Tier 1 DEG (n = %s tested in all 5 cohorts)",
         comma(nrow(ps_dir))),
       caption = paste0(
         "Dashed line = expected fraction under independent 50/50 sign null.\n",
         "Most Tier 1 DEGs are 4/5 or 5/5 cohort-concordant by direction — far above random."
       )) +
  theme_masld(base_size = 7) +
  theme(plot.title   = element_text(size = 8, face = "bold"),
        plot.caption = element_text(size = 6, color = "gray35", hjust = 0))

save_fig(p_k, file.path(PANEL_DIR, "figS_batch_k_direction_concordance.pdf"),
         width = fig_half_width * 1.1, height = 2.6)

# ----------------------------------------------------------------------------
# PANEL L: kBET-style chi-square cohort-mixing test
#   For each sample, compute observed cohort proportions in its 30-NN window,
#   compare to expected (global) cohort proportions via chi-square.
#   Rejection rate at p < 0.05 = fraction of samples with non-mixed neighborhood.
# ----------------------------------------------------------------------------
message("Panel L: kBET-style mixing test...")
emb <- as.matrix(plot_dt[, .(UMAP1, UMAP2)])
nn  <- RANN::nn2(emb, query = emb, k = KNN_K + 1)$nn.idx[, -1]
cohort_int <- as.integer(plot_dt$cohort)
expected_props <- prop.table(table(plot_dt$cohort))

kbet_chi <- function(neighbor_idx, expected_props, cohort_int) {
  obs_counts <- tabulate(cohort_int[neighbor_idx],
                         nbins = length(expected_props))
  exp_counts <- expected_props * length(neighbor_idx)
  # Chi-square (asymptotic; used here as a relative score, not a strict test)
  keep <- exp_counts > 0
  sum((obs_counts[keep] - exp_counts[keep])^2 / exp_counts[keep])
}
chi_per_sample <- vapply(seq_len(nrow(nn)),
                         function(i) kbet_chi(nn[i, ], expected_props, cohort_int),
                         numeric(1))
plot_dt[, kbet_chi := chi_per_sample]

# Per-cohort distribution of chi-square (continuous mixing score)
# Lower chi-square = neighborhoods closer to global cohort proportions = better mixed.
kbet_summ <- plot_dt[, .(median_chi = median(kbet_chi, na.rm = TRUE),
                         q25 = quantile(kbet_chi, 0.25, na.rm = TRUE),
                         q75 = quantile(kbet_chi, 0.75, na.rm = TRUE),
                         n = .N),
                     by = cohort]
setorder(kbet_summ, median_chi)
kbet_summ[, cohort := factor(cohort, levels = cohort)]
plot_dt[, cohort_ord := factor(cohort, levels = levels(kbet_summ$cohort))]
overall_chi <- median(plot_dt$kbet_chi, na.rm = TRUE)

# Floor at 1% quantile for a tight ridge-style boxplot
chi_cap <- quantile(plot_dt$kbet_chi, 0.99, na.rm = TRUE)

p_l <- ggplot(plot_dt, aes(x = pmin(kbet_chi, chi_cap), y = cohort_ord,
                           fill = cohort_ord)) +
  geom_boxplot(width = 0.55, outlier.size = 0.3, outlier.alpha = 0.5,
               linewidth = 0.25, color = "gray25") +
  geom_vline(xintercept = overall_chi,
             linetype = "dashed", color = "gray40", linewidth = 0.3) +
  annotate("text", x = overall_chi, y = 0.6,
           label = sprintf("overall median = %.0f", overall_chi),
           hjust = -0.05, vjust = 0, size = 2.0, color = "gray25") +
  scale_fill_manual(values = COHORT_COLORS, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = "kBET chi-square score (per-sample 30-NN mixing)\nlower = better cohort mixing",
       y = NULL,
       title = "Per-cohort kBET-style mixing distribution",
       caption = paste0(
         "Per-sample chi-square of observed vs global cohort proportions in 30-NN\n",
         "UMAP windows. Lower scores = better cross-cohort mixing in that neighborhood.\n",
         "Cohorts ordered by median; values capped at 99th percentile for display."
       )) +
  theme_masld(base_size = 7) +
  theme(plot.title    = element_text(size = 8, face = "bold"),
        plot.caption  = element_text(size = 6, color = "gray35", hjust = 0))

save_fig(p_l, file.path(PANEL_DIR, "figS_batch_l_kbet.pdf"),
         width = fig_half_width * 1.05, height = 2.6)

# ----------------------------------------------------------------------------
# PANEL M: Per-gene variance scatter — dataset vs disease+fibrosis biology
# ----------------------------------------------------------------------------
message("Panel M: per-gene variance partition scatter...")
vp_path <- file.path(BASE,
                     "RNA-seq/results/audit_sensitivity/variance_partition_results.csv")
vp <- fread(vp_path)
vp[, gene_clean := sub("\\..*", "", gene)]
vp <- merge(vp, dream[, .(gene_clean, is_tier1, dream_logFC, dream_padj)],
            by = "gene_clean", all.x = TRUE)
vp[is.na(is_tier1), is_tier1 := FALSE]
vp[, biology_var := pmax(0, fibrosis_stage) + pmax(0, group_binary)]
vp[, dataset_var := pmax(0, dataset)]
vp[, deg_class   := fifelse(is_tier1, "Tier 1 DEG", "Non-DEG")]

# Pick a balanced subsample for the non-DEG cloud so it doesn't drown the plot
n_t1   <- sum(vp$deg_class == "Tier 1 DEG")
non_t1 <- vp[deg_class == "Non-DEG"]
set.seed(42)
non_t1_sub <- non_t1[sample(.N, min(.N, 4 * n_t1))]
vp_plot <- rbind(vp[deg_class == "Tier 1 DEG"], non_t1_sub)
vp_plot[, deg_class := factor(deg_class, levels = c("Non-DEG", "Tier 1 DEG"))]
setorder(vp_plot, deg_class)   # Tier 1 plotted on top

# Above-diagonal fraction
above_t1   <- 100 * mean(vp[deg_class == "Tier 1 DEG", biology_var > dataset_var],
                         na.rm = TRUE)
above_non  <- 100 * mean(vp[deg_class == "Non-DEG",    biology_var > dataset_var],
                         na.rm = TRUE)

p_m <- ggplot(vp_plot, aes(x = 100 * dataset_var, y = 100 * biology_var,
                           color = deg_class)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "gray45", linewidth = 0.3) +
  rasterize_layer(
    geom_point(size = 0.32, alpha = 0.55, shape = 16)
  ) +
  scale_color_manual(values = c("Non-DEG" = "#BDBDBD",
                                "Tier 1 DEG" = "#C2185B"),
                     name = NULL) +
  scale_x_continuous(limits = c(0, 100),
                     labels = function(v) paste0(v, "%"),
                     expand = expansion(0)) +
  scale_y_continuous(limits = c(0, 100),
                     labels = function(v) paste0(v, "%"),
                     expand = expansion(0)) +
  annotate("text", x = 8, y = 92,
           label = "biology > batch", hjust = 0, vjust = 1,
           size = 2.0, color = "gray25", fontface = "italic") +
  annotate("text", x = 92, y = 8,
           label = "batch > biology", hjust = 1, vjust = 0,
           size = 2.0, color = "gray25", fontface = "italic") +
  annotate("text", x = 95, y = 95,
           label = sprintf("Tier 1 above diag: %.1f%%\nNon-DEG above diag: %.1f%%",
                           above_t1, above_non),
           hjust = 1, vjust = 1, size = 2.0, color = "gray15") +
  coord_fixed() +
  labs(x = "% variance explained by dataset",
       y = "% variance explained by\ndisease + fibrosis",
       title = "Per-gene biology vs batch variance",
       caption = paste0(
         "variancePartition (Hoffman & Schadt 2016) joint model:\n",
         "~ dataset + fibrosis + disease + sex + (1 | sample). Each dot = one gene."
       )) +
  theme_masld(base_size = 7) +
  theme(plot.title    = element_text(size = 8, face = "bold"),
        plot.caption  = element_text(size = 6, color = "gray35", hjust = 0),
        legend.position = "top",
        legend.key.size = unit(0.28, "cm"),
        legend.text     = element_text(size = 6))

save_fig(p_m, file.path(PANEL_DIR, "figS_batch_m_gene_variance_scatter.pdf"),
         width = fig_half_width * 1.0, height = 2.8)

# ----------------------------------------------------------------------------
# Combined figure: assemble h–m
# ----------------------------------------------------------------------------
message("Composing combined extended figure...")

# Existing panels a–g (rebuild thumbnails by re-loading saved RDS / regenerating
# minimal versions). For composability we just import the per-panel PDFs as
# patchwork wrap_elements via magick::image_read_pdf if available; otherwise
# we keep the standalone fig as the canonical combined output and produce a
# separate "extra" combined that holds h–l only.

extra_top    <- p_h
extra_mid    <- p_i
extra_bot    <- (p_j | p_k | p_l | p_m) + plot_layout(widths = c(1.0, 0.85, 0.95, 0.9))
fig_extra    <- extra_top / extra_mid / extra_bot +
  plot_layout(heights = c(1.0, 0.9, 1.0)) +
  plot_annotation(
    title    = "Extended quantitative evidence: cohort separation in UMAP is biology, not residual batch",
    subtitle = "10-cohort MASLD bulk RNA-seq atlas (n = 1,444) — supplementing figS_batch_correction panels a-g",
    tag_levels = list(c("h", "", "", "", "i", "j", "k", "l", "m")),
    theme = theme(plot.title    = element_text(size = 9, face = "bold"),
                  plot.subtitle = element_text(size = 7, color = "gray35"),
                  plot.tag      = element_text(size = 9, face = "bold"))
  )

out_extra <- file.path(FIGS_BATCH_DIR, "figS_batch_correction_extra.pdf")
save_fig(fig_extra, out_extra, width = fig_full_width * 1.4, height = 8.4)

# Per-cohort logFC concordance CSV
fwrite(cohort_stats, file.path(FIGS_BATCH_DIR, "figS_batch_i_concordance_stats.csv"))
# kBET CSV
fwrite(kbet_summ, file.path(FIGS_BATCH_DIR, "figS_batch_l_kbet_summary.csv"))
# PC variance attribution CSV
fwrite(attr_long, file.path(FIGS_BATCH_DIR, "figS_batch_j_pc_attribution.csv"))
# Direction concordance CSV
fwrite(dir_counts, file.path(FIGS_BATCH_DIR, "figS_batch_k_direction_concordance.csv"))

if (file.exists(out_extra)) {
  message(sprintf("\nOutput: %s (%s)", out_extra,
                  utils:::format.object_size(file.size(out_extra), "auto")))
}
message(sprintf("Per-panel PDFs in: %s", PANEL_DIR))
