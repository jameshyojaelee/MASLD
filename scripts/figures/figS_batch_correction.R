#!/usr/bin/env Rscript
# figS_batch_correction.R — Six-panel supplementary figure addressing whether
# the cohort-correlated structure visible in the integrated UMAP reflects
# residual batch effect or biology.
#
# Headline argument:
#   The dream mega-analysis uses (1 | dataset) as a random intercept; cohort
#   variance lives in the random effect rather than being scrubbed from the
#   expression matrix. UMAP, being a low-rank visualization, surfaces residual
#   technical AND biological cohort structure equally. They are different.
#
# Panels (each saved to panels/ as well):
#   A: Cohort × disease state UMAP pair (existing UMAP coords, paired layout)
#   B: Cohort × diagnosis confounding heatmap (counts + %)
#   C: kNN-based mixing — iLISI-style metric for cohort vs disease label
#   D: Variance partition — fraction of per-gene variance attributable to
#      dataset / disease / fibrosis / sex / residual, stratified by Tier 1
#      DEG vs non-DEG
#   E: LOO-CV recovery of integrated DEGs by held-out cohort (downstream
#      stability check; uses fig1g LOO data at |LFC|=0.5)
#   F: "UMAP ≠ DE" — schematic / text panel explaining the distinction
#
# Output:
#   figures/supplementary/figS_methods_validation/batch_correction/figS_batch_correction.pdf
#   figures/supplementary/figS_methods_validation/batch_correction/panels/figS_batch_{a..f}.pdf

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(RANN)
})

set.seed(42)

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR <- file.path(FIGS_BATCH_DIR, "panels")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)

# Match fig1c first-author labels
COHORT_LABEL <- c(
  GSE126848   = "Suppli",     GSE130970   = "Hoang",
  GSE135251   = "Govaere",    GSE162694   = "Bril",
  GSE167523   = "Kozumi",     GSE174478   = "Kawamura",
  GSE193066   = "Hoshida",    GSE213621   = "Chen",
  GSE240729   = "Verschuren", PRJNA512027 = "Gerhard"
)
COHORT_COLORS <- c(
  "Suppli"="#1F77B4","Hoang"="#FF7F0E","Govaere"="#2CA02C","Bril"="#D62728",
  "Kozumi"="#9467BD","Kawamura"="#8C564B","Hoshida"="#E377C2","Chen"="#7F7F7F",
  "Verschuren"="#BCBD22","Gerhard"="#17BECF"
)

PADJ_INT <- 0.05
LFC_INT  <- 0.5    # Tier 1
KNN_K    <- 30      # neighborhood size for the mixing metric

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
message(sprintf("Tier 1 DEGs (padj<%.2g, |LFC|>%.1f): %s",
                PADJ_INT, LFC_INT, comma(n_tier1)))

# ----------------------------------------------------------------------------
# PANEL A: Cohort + disease UMAP pair
# ----------------------------------------------------------------------------
message("Building panel A (UMAP cohort/disease pair)...")
disease_colors <- c(Control = masld_colors$control,
                    Disease = masld_colors$nash)

base_umap <- function(dt, color_var, palette, title, legend_ncol = 2) {
  ggplot(dt, aes(x = UMAP1, y = UMAP2, color = .data[[color_var]])) +
    rasterize_layer(geom_point(size = 0.32, alpha = 0.7, shape = 16)) +
    scale_color_manual(values = palette, name = NULL,
                       na.value = "gray85", drop = FALSE) +
    coord_fixed() +
    guides(color = guide_legend(ncol = legend_ncol,
                                override.aes = list(size = 1.2, alpha = 1))) +
    labs(x = "UMAP 1", y = "UMAP 2", title = title) +
    theme_masld(base_size = 7) +
    theme(plot.title = element_text(size = 8, face = "bold"),
          legend.text = element_text(size = 6),
          legend.key.size = unit(0.22, "cm"),
          legend.position = "right")
}

p_a_cohort  <- base_umap(plot_dt, "cohort", COHORT_COLORS, "by cohort", legend_ncol = 2)
p_a_disease <- base_umap(plot_dt, "disease_state", disease_colors, "by disease state", legend_ncol = 1)
p_a <- p_a_cohort | p_a_disease

save_fig(p_a, file.path(PANEL_DIR, "figS_batch_a_umap_pair.pdf"),
         width = fig_full_width * 0.85, height = 2.6)

# ----------------------------------------------------------------------------
# PANEL B: Cohort × disease confounding heatmap
# ----------------------------------------------------------------------------
message("Building panel B (cohort × diagnosis confounding)...")
plot_dt[, dx := fcase(
  diagnosis_harmonized == "Control",                                  "Control",
  diagnosis_harmonized == "NAFL",                                     "NAFL/MASL",
  diagnosis_harmonized %in% c("Borderline", "NASH"),                  "NASH/MASH",
  default = "Unstaged")]
plot_dt[, dx := factor(dx,
                       levels = c("Control", "NAFL/MASL", "NASH/MASH", "Unstaged"))]

confound <- plot_dt[, .N, by = .(cohort, dx)]
confound[, total := sum(N), by = cohort]
confound[, pct  := 100 * N / total]
confound[, label := ifelse(N > 0, sprintf("%d\n(%.0f%%)", N, pct), "")]
# Order cohorts by total
order_dt <- plot_dt[, .N, by = cohort][order(-N)]
confound[, cohort := factor(cohort, levels = rev(order_dt$cohort))]

p_b <- ggplot(confound, aes(x = dx, y = cohort, fill = pct)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = label), size = 1.95, color = "gray15", lineheight = 0.85) +
  scale_fill_gradientn(
    colors = c("#FFFFFF", "#FCE4EC", "#F48FB1", "#C2185B"),
    limits = c(0, 100), name = "%",
    labels = function(v) paste0(v, "%")) +
  labs(x = NULL, y = NULL,
       title = "Cohort × disease-state confounding (sample counts)") +
  theme_masld(base_size = 7) +
  theme(plot.title = element_text(size = 8, face = "bold"),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.key.height = unit(0.32, "cm"),
        legend.key.width  = unit(0.18, "cm"),
        legend.position = "right")

save_fig(p_b, file.path(PANEL_DIR, "figS_batch_b_confounding.pdf"),
         width = fig_half_width * 1.0, height = 2.6)

# ----------------------------------------------------------------------------
# PANEL C: kNN-based mixing metric (iLISI-style)
#   For each sample, compute the inverse Simpson's index of label proportions
#   among its KNN_K nearest neighbors in UMAP space.
#   - cohort iLISI: high = better cross-cohort mixing (max = n_cohorts = 10)
#   - disease iLISI: lower means cleaner separation between Control/Disease
#     (max = 2, so values close to 2 = mixed = bad biology preservation)
# ----------------------------------------------------------------------------
message("Building panel C (kNN mixing metric)...")
ilisi <- function(coords, labels, k) {
  nn  <- RANN::nn2(coords, query = coords, k = k + 1)$nn.idx[, -1]
  lab <- as.integer(as.factor(labels))
  vapply(seq_len(nrow(nn)), function(i) {
    p <- table(lab[nn[i, ]]) / k
    1 / sum(p^2)
  }, numeric(1))
}

emb <- as.matrix(plot_dt[, .(UMAP1, UMAP2)])
plot_dt[, lisi_cohort  := ilisi(emb, cohort,         KNN_K)]
plot_dt[, lisi_disease := ilisi(emb, disease_state,  KNN_K)]

n_coh <- length(unique(plot_dt$cohort))
n_dis <- length(unique(plot_dt$disease_state))

lisi_long <- rbind(
  plot_dt[, .(metric = "Cohort iLISI",   value = lisi_cohort,   max = n_coh)],
  plot_dt[, .(metric = "Disease iLISI",  value = lisi_disease,  max = n_dis)]
)
lisi_long[, metric := factor(metric, levels = c("Cohort iLISI", "Disease iLISI"))]

# Per-metric medians for annotation
lisi_summary <- lisi_long[, .(med = median(value, na.rm = TRUE),
                              q25 = quantile(value, 0.25, na.rm = TRUE),
                              q75 = quantile(value, 0.75, na.rm = TRUE),
                              max = unique(max)), by = metric]

p_c <- ggplot(lisi_long, aes(x = metric, y = value, fill = metric)) +
  geom_violin(trim = FALSE, color = NA, alpha = 0.7, scale = "width", width = 0.85) +
  geom_boxplot(width = 0.12, fill = "white", color = "gray25",
               outlier.shape = NA, linewidth = 0.3) +
  geom_segment(data = lisi_summary,
               aes(x = as.numeric(metric) - 0.45, xend = as.numeric(metric) + 0.45,
                   y = max, yend = max),
               linetype = "dashed", color = "gray40", linewidth = 0.3,
               inherit.aes = FALSE) +
  geom_text(data = lisi_summary,
            aes(x = metric, y = max + 0.18,
                label = sprintf("max possible\n(%d unique labels)", max)),
            inherit.aes = FALSE, size = 1.7, color = "gray35", lineheight = 0.85) +
  geom_text(data = lisi_summary,
            aes(x = metric, y = q25 - 0.4,
                label = sprintf("median = %.2f", med)),
            inherit.aes = FALSE, size = 2.1, fontface = "bold", color = "gray15") +
  scale_fill_manual(values = c("Cohort iLISI" = "#9E9E9E",
                               "Disease iLISI" = "#C2185B"),
                    guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.08))) +
  labs(x = NULL, y = "iLISI",
       title = sprintf("Local label diversity (k = %d nearest neighbors)", KNN_K),
       caption = paste0(
         "Higher Cohort iLISI = better cohort mixing in the embedding.\n",
         "Lower Disease iLISI = cleaner Control/Disease separation (biology preserved)."
       )) +
  theme_masld(base_size = 7) +
  theme(plot.title    = element_text(size = 8, face = "bold"),
        plot.caption  = element_text(size = 6, color = "gray35", hjust = 0))

save_fig(p_c, file.path(PANEL_DIR, "figS_batch_c_ilisi.pdf"),
         width = fig_half_width * 0.9, height = 2.6)

# ----------------------------------------------------------------------------
# PANEL D: Variance partition by DEG class
# ----------------------------------------------------------------------------
message("Building panel D (variance partition)...")
vp_path <- file.path(BASE, "RNA-seq/results/audit_sensitivity/variance_partition_results.csv")
vp <- fread(vp_path)
vp[, gene_clean := sub("\\..*", "", gene)]
vp <- merge(vp,
            dream[, .(gene_clean, is_tier1)],
            by = "gene_clean", all.x = TRUE)
vp[is.na(is_tier1), is_tier1 := FALSE]
vp[, deg_class := fifelse(is_tier1, "Tier 1 DEG", "Non-DEG")]

# Long form
vp_long <- melt(vp,
                id.vars = c("gene_clean", "deg_class"),
                measure.vars = c("dataset", "group_binary",
                                 "fibrosis_stage", "sex_covar", "Residuals"),
                variable.name = "covariate", value.name = "frac")
vp_long[, frac := pmin(pmax(frac, 0), 1)]
covar_label <- c(
  dataset = "Dataset",
  group_binary = "Disease (binary)",
  fibrosis_stage = "Fibrosis stage",
  sex_covar = "Sex",
  Residuals = "Residual"
)
vp_long[, covariate := factor(covar_label[as.character(covariate)],
                              levels = covar_label)]

vp_summ <- vp_long[, .(median_pct = 100 * median(frac, na.rm = TRUE)),
                   by = .(covariate, deg_class)]

p_d <- ggplot(vp_long, aes(x = covariate, y = 100 * frac, fill = deg_class)) +
  geom_violin(position = position_dodge(width = 0.75),
              width = 0.7, scale = "width", color = NA, alpha = 0.85,
              trim = FALSE) +
  geom_boxplot(position = position_dodge(width = 0.75),
               width = 0.18, fill = "white", color = "gray25",
               outlier.shape = NA, linewidth = 0.25) +
  scale_fill_manual(values = c("Non-DEG" = "#9E9E9E",
                               "Tier 1 DEG" = "#C2185B"),
                    name = NULL) +
  scale_y_continuous(labels = function(v) paste0(v, "%"),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "% variance per gene",
       title = "Per-gene variance partition (Tier 1 DEG vs non-DEG)",
       caption = paste0(
         "Dataset variance is high in *non-DEGs* (median = ",
         round(vp_summ[covariate == "Dataset" & deg_class == "Non-DEG", median_pct], 1),
         "%) but biology dominates Tier 1 DEGs."
       )) +
  theme_masld(base_size = 7) +
  theme(plot.title   = element_text(size = 8, face = "bold"),
        plot.caption = element_text(size = 6, color = "gray35", hjust = 0),
        legend.position = "top",
        legend.key.size = unit(0.28, "cm"))

save_fig(p_d, file.path(PANEL_DIR, "figS_batch_d_varpart.pdf"),
         width = fig_full_width * 0.55, height = 2.6)

# ----------------------------------------------------------------------------
# PANEL E: LOO-CV stability of integrated DEGs at |LFC|=0.5
# ----------------------------------------------------------------------------
message("Building panel E (LOO-CV recovery)...")
loo_path <- file.path(FIG1_DIR, "panels", "fig1g_loo_stability_data.csv")
loo <- fread(loo_path)
loo_05 <- loo[lfc_cutoff == 0.5]

# Apply first-author labels
loo_05[, cohort_label := COHORT_LABEL[held_out]]
setorder(loo_05, pct_recovered)
loo_05[, cohort_label := factor(cohort_label, levels = cohort_label)]

mean_rec <- mean(loo_05$pct_recovered)
mean_jac <- mean(loo_05$jaccard)

p_e <- ggplot(loo_05, aes(x = pct_recovered, y = cohort_label)) +
  geom_segment(aes(x = 0, xend = pct_recovered, yend = cohort_label),
               color = "gray70", linewidth = 0.4) +
  geom_point(aes(size = pct_held), color = "#C2185B", alpha = 0.85) +
  geom_text(aes(label = sprintf("%.0f%%", pct_recovered)),
            hjust = -0.3, size = 1.95, color = "gray20") +
  geom_vline(xintercept = mean_rec, linetype = "dashed",
             color = "gray40", linewidth = 0.3) +
  annotate("text",
           x = mean_rec, y = 0.6,
           label = sprintf("mean recovery = %.1f%%\nJaccard = %.2f",
                           mean_rec, mean_jac),
           hjust = -0.05, vjust = 0, size = 2.1, color = "gray25",
           lineheight = 0.85) +
  scale_size_continuous(range = c(1.2, 3.4),
                        breaks = c(5, 15, 30),
                        name = "% patients\nheld out") +
  scale_x_continuous(limits = c(0, 110),
                     breaks = c(0, 25, 50, 75, 100),
                     labels = function(v) paste0(v, "%"),
                     expand = expansion(mult = 0)) +
  labs(x = "% Tier 1 DEGs recovered when cohort held out",
       y = NULL,
       title = "LOO-CV recovery of integrated DEGs (|LFC| > 0.5)") +
  theme_masld(base_size = 7) +
  theme(plot.title = element_text(size = 8, face = "bold"),
        legend.position  = "right",
        legend.key.size  = unit(0.3, "cm"),
        legend.title     = element_text(size = 6.5),
        legend.text      = element_text(size = 6))

save_fig(p_e, file.path(PANEL_DIR, "figS_batch_e_loo_recovery.pdf"),
         width = fig_half_width * 1.05, height = 2.6)

# ----------------------------------------------------------------------------
# PANEL F: "UMAP ≠ DE" schematic / text panel
# ----------------------------------------------------------------------------
message("Building panel F (UMAP vs DE schematic)...")
text_lines <- list(
  list(label = "UMAP",
       y = 0.92, size = 5.0, face = "bold", color = "#1565C0"),
  list(label = "Low-rank 2D embedding for visualization.",
       y = 0.84, size = 2.6, face = "plain", color = "gray25"),
  list(label = "Operates on a denoised PCA representation of",
       y = 0.79, size = 2.6, face = "plain", color = "gray25"),
  list(label = "the full expression matrix; surfaces residual",
       y = 0.74, size = 2.6, face = "plain", color = "gray25"),
  list(label = "technical AND biological cohort structure.",
       y = 0.69, size = 2.6, face = "plain", color = "gray25"),
  list(label = "Cohort × disease confounding looks like batch.",
       y = 0.62, size = 2.6, face = "italic", color = "#1565C0"),

  list(label = "dream mega-analysis",
       y = 0.45, size = 5.0, face = "bold", color = "#C2185B"),
  list(label = "Mixed-effects model on the full gene matrix:",
       y = 0.37, size = 2.6, face = "plain", color = "gray25"),
  list(label = "  ~ disease + sex + (1 | dataset)",
       y = 0.31, size = 2.6, face = "plain", color = "gray15"),
  list(label = "Cohort variance lives in a random intercept;",
       y = 0.25, size = 2.6, face = "plain", color = "gray25"),
  list(label = "fixed effects extract the disease signal across",
       y = 0.20, size = 2.6, face = "plain", color = "gray25"),
  list(label = "cohorts without scrubbing it from expression.",
       y = 0.15, size = 2.6, face = "plain", color = "gray25"),
  list(label = "DE call does not depend on UMAP geometry.",
       y = 0.08, size = 2.6, face = "italic", color = "#C2185B")
)
text_dt <- rbindlist(lapply(text_lines, as.data.table))

p_f <- ggplot(text_dt, aes(x = 0.5, y = y, label = label)) +
  geom_text(aes(size = size, fontface = face, color = color),
            hjust = 0.5, lineheight = 0.85) +
  scale_size_identity() +
  scale_color_identity() +
  scale_x_continuous(limits = c(0, 1), expand = expansion(0)) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(0)) +
  labs(title = "UMAP visualization ≠ DE inference") +
  theme_void(base_size = 7) +
  theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5,
                                  color = "gray10"),
        plot.background = element_rect(fill = "#FAFAFA", color = "gray85",
                                       linewidth = 0.3),
        plot.margin = margin(8, 8, 6, 8))

save_fig(p_f, file.path(PANEL_DIR, "figS_batch_f_umap_vs_de.pdf"),
         width = fig_half_width * 0.95, height = 2.6)

# ----------------------------------------------------------------------------
# Combined figure
# ----------------------------------------------------------------------------
message("Composing combined figure...")
top_row    <- p_a
mid_row    <- p_b | p_c | p_d
bot_row    <- p_e | p_f

fig <- top_row / mid_row / bot_row +
  plot_layout(heights = c(1, 1, 1)) +
  plot_annotation(
    title    = "Batch correction adequacy: cohort structure in UMAP is biology, not residual batch",
    subtitle = "10 cohorts | 1,444 QC-passing samples | dream uses (1 | dataset) random intercept",
    tag_levels = "a",
    theme = theme(plot.title    = element_text(size = 9, face = "bold"),
                  plot.subtitle = element_text(size = 7, color = "gray35"),
                  plot.tag      = element_text(size = 9, face = "bold"))
  )

out_pdf <- file.path(FIGS_BATCH_DIR, "figS_batch_correction.pdf")
save_fig(fig, out_pdf, width = fig_full_width * 1.4, height = 8.4)

# Per-sample iLISI export
fwrite(plot_dt[, .(sample_id, dataset, cohort, disease_state,
                   UMAP1, UMAP2, lisi_cohort, lisi_disease)],
       file.path(FIGS_BATCH_DIR, "figS_batch_data.csv"))

if (file.exists(out_pdf)) {
  message(sprintf("\nOutput: %s (%s)", out_pdf,
                  utils:::format.object_size(file.size(out_pdf), "auto")))
  message(sprintf("Per-panel PDFs: %s/figS_batch_{a..f}_*.pdf", PANEL_DIR))
}
