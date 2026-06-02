#!/usr/bin/env Rscript
# figS_pi_concordance_panels.R
# Four additional PI-facing panels proving that the population signature is
# recapitulated at the individual level: healthy controls → 0, disease → up.
#
# Output — same patient_concordance/ subfolder as figS_patient_concordance.R:
#   E1_fibrosis_stage_score.pdf   — signature score by fibrosis stage (ctrl → F4)
#   E2_waterfall.pdf              — all samples ranked by signature score
#   E3_roc_curve.pdf              — ROC curve: disease vs control discrimination
#   E4_cohort_boxplots.pdf        — per-cohort Control vs Disease violin + box
#
# Canonical cutoff: dream |LFC| > 0.5 (padj < 0.05).
# Data: same sources as figS_patient_concordance.R (no new inputs).

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(scales)
  library(edgeR)
  library(yaml)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))

OUTDIR      <- file.path(BASE, "figures/supplementary/figS_lfc_sensitivity/patient_concordance")
INT_RESULTS <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")

# ── Publication theme (matches existing panels) ───────────────────────────────
theme_pub <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line        = element_line(colour = "black", linewidth = 0.3),
    axis.ticks       = element_line(colour = "black", linewidth = 0.3),
    legend.background = element_blank(),
    legend.key        = element_blank(),
    strip.background  = element_blank(),
    strip.text        = element_text(face = "bold", size = 10),
    plot.title        = element_text(face = "bold", size = 12),
    plot.subtitle     = element_text(size = 9, color = "grey40"),
    axis.title        = element_text(size = 10),
    axis.text         = element_text(size = 9),
    legend.text       = element_text(size = 9),
    legend.title      = element_text(size = 9),
    plot.margin       = margin(8, 10, 8, 8)
  )
theme_set(theme_pub)

# ── Colors ────────────────────────────────────────────────────────────────────
col_disease <- "#00695C"   # teal (matches existing panels in patient_concordance/)
col_control <- "#9E9E9E"   # neutral gray (project rule: controls must be gray)
# fibrosis_stage_colors from publication_theme.R:
#   F0="#E3F2FD"  F1="#90CAF9"  F2="#42A5F5"  F3="#1565C0"  F4="#0D47A1"
grp_colors <- c("Disease" = col_disease, "Control" = col_control)

# ════════════════════════════════════════════════════════════════════════════
# Data loading
# ════════════════════════════════════════════════════════════════════════════

message("Loading dream STAR results...")
dream <- fread(file.path(INT_RESULTS, "dream_results_ashr.csv"))
setnames(dream, "logFC", "dream_logFC", skip_absent = FALSE)
setnames(dream, "padj",  "dream_padj",  skip_absent = FALSE)

message("Loading patient LFC matrix (disease samples)...")
lfc_mat    <- fread(file.path(INT_RESULTS, "patient_lfc_matrix.csv.gz"))
setkey(lfc_mat, gene)
sample_ids <- setdiff(colnames(lfc_mat), "gene")
message(sprintf("  %d genes x %d disease patients", nrow(lfc_mat), length(sample_ids)))

message("Loading unified metadata...")
meta <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv"))

message("Loading merged DGE for control LFC computation...")
dge_s  <- readRDS(file.path(INT_RESULTS, "merged_dge.rds"))
ycfg_d <- yaml::read_yaml(file.path(BASE, "config/human_datasets.yaml"))$datasets
mega_c <- names(Filter(function(d) isTRUE(d$de$include_in_mega), ycfg_d))
dge_s  <- dge_s[, dge_s$samples$dataset %in% mega_c]
lcpm_s <- edgeR::cpm(dge_s, log = TRUE, prior.count = 1)

sinfo    <- data.table(sample_id = colnames(dge_s),
                       dataset   = dge_s$samples$dataset,
                       group     = as.character(dge_s$samples$group_binary))
valid_ds <- sinfo[, .(nc = sum(group == "Control"), nd = sum(group == "Disease")),
                  by = dataset][nc > 0 & nd > 0, dataset]
message(sprintf("  Canonical cohorts with controls: %s", paste(valid_ds, collapse = ", ")))

# Within-dataset control means (same reference as disease samples)
cmeans <- lapply(setNames(valid_ds, valid_ds), function(ds) {
  ids <- sinfo[dataset == ds & group == "Control", sample_id]
  if (length(ids) == 1) lcpm_s[, ids] else rowMeans(lcpm_s[, ids, drop = FALSE])
})

ctrl_sids <- sinfo[group == "Control" & dataset %in% valid_ds, sample_id]
ctrl_lfc  <- matrix(NA_real_, nrow = nrow(lcpm_s), ncol = length(ctrl_sids),
                    dimnames = list(rownames(lcpm_s), ctrl_sids))
for (ds in valid_ds) {
  ids <- sinfo[dataset == ds & group == "Control", sample_id]
  for (sid in ids) ctrl_lfc[, sid] <- lcpm_s[, sid] - cmeans[[ds]]
}
message(sprintf("  Control LFC matrix: %d genes x %d controls",
                nrow(ctrl_lfc), ncol(ctrl_lfc)))

# ── Align all objects to common gene set ─────────────────────────────────────
common_genes <- Reduce(intersect, list(dream$gene, lfc_mat$gene, rownames(ctrl_lfc)))
message(sprintf("  Common genes: %d", length(common_genes)))

dream_sub <- dream[gene %in% common_genes]
lfc_sub   <- lfc_mat[gene %in% common_genes]
setkey(lfc_sub, gene)
lfc_sub   <- lfc_sub[dream_sub$gene]  # row order matches dream_sub
lfc_m     <- as.matrix(lfc_sub[, sample_ids, with = FALSE])
rownames(lfc_m) <- dream_sub$gene
ctrl_lfc  <- ctrl_lfc[dream_sub$gene, ]

# ── Signature score at canonical cutoff ──────────────────────────────────────
CUT <- 0.5
up_g    <- dream_sub$gene[dream_sub$dream_padj < 0.05 & dream_sub$dream_logFC >  CUT]
dn_g    <- dream_sub$gene[dream_sub$dream_padj < 0.05 & dream_sub$dream_logFC < -CUT]
n_total <- length(up_g) + length(dn_g)
message(sprintf("  DEGs at |LFC|>%.1f: %d up + %d down = %d total",
                CUT, length(up_g), length(dn_g), n_total))

sig_score <- function(mat, up, dn, n) {
  su  <- if (length(up)) colSums(mat[up, , drop = FALSE]) else rep(0, ncol(mat))
  sd_ <- if (length(dn)) colSums(mat[dn, , drop = FALSE]) else rep(0, ncol(mat))
  (su - sd_) / n
}

dis_scores  <- sig_score(lfc_m,    up_g, dn_g, n_total)
ctrl_scores <- sig_score(ctrl_lfc, up_g, dn_g, n_total)

# ── Per-sample score table with metadata ─────────────────────────────────────
scores_dis <- merge(
  data.table(sample_id = names(dis_scores), sig_score = as.numeric(dis_scores),
             group_binary = "Disease"),
  meta[, .(sample_id, dataset, fibrosis_stage = as.character(fibrosis_stage))],
  by = "sample_id", all.x = TRUE)

scores_ctrl <- merge(
  data.table(sample_id = names(ctrl_scores), sig_score = as.numeric(ctrl_scores),
             group_binary = "Control"),
  meta[, .(sample_id, dataset, fibrosis_stage = as.character(fibrosis_stage))],
  by = "sample_id", all.x = TRUE)

scores_all <- rbind(scores_dis, scores_ctrl)

message(sprintf("\nAt canonical |LFC|>%.1f:", CUT))
message(sprintf("  Disease  n=%d  median=%.3f", nrow(scores_dis),
                median(scores_dis$sig_score, na.rm = TRUE)))
message(sprintf("  Control  n=%d  median=%.3f", nrow(scores_ctrl),
                median(scores_ctrl$sig_score, na.rm = TRUE)))

# ════════════════════════════════════════════════════════════════════════════
# E1: Mean log₂FC in disease direction by fibrosis stage
# Metric: for each patient, mean( patient_log2FC ) over all DEGs after
# flipping the sign of down-regulated genes so that "up = disease direction".
# Equivalent to the signature score; values are in log₂FC units.
# Controls = the reference mean → LFC ≈ 0.  Disease rises F0→F4.
# ════════════════════════════════════════════════════════════════════════════
message("\n── E1: Mean log2FC (disease direction) by fibrosis stage ──")

mean_lfc_dis  <- as.numeric(sig_score(lfc_m,    up_g, dn_g, n_total))
mean_lfc_ctrl <- as.numeric(sig_score(ctrl_lfc, up_g, dn_g, n_total))
message(sprintf("  Mean log2FC disease direction — disease median=%.3f, control median=%.3f",
                median(mean_lfc_dis), median(mean_lfc_ctrl)))

scores_dis[,  mean_lfc_up := mean_lfc_dis]
scores_ctrl[, mean_lfc_up := mean_lfc_ctrl]
scores_all <- rbind(scores_dis, scores_ctrl)

# Build stage groups (missing stage excluded from E1 only)
scores_stage <- copy(scores_all)
scores_stage[, stage_group := ifelse(
  group_binary == "Control",
  "Control",
  ifelse(fibrosis_stage %in% as.character(0:4),
         paste0("F", fibrosis_stage),
         NA_character_)
)]
scores_stage <- scores_stage[!is.na(stage_group)]
scores_stage[, stage_group := factor(stage_group,
                                     levels = c("Control", paste0("F", 0:4)))]

stage_fill <- c("Control" = col_control,
                "F0"      = unname(fibrosis_stage_colors["F0"]),
                "F1"      = unname(fibrosis_stage_colors["F1"]),
                "F2"      = unname(fibrosis_stage_colors["F2"]),
                "F3"      = unname(fibrosis_stage_colors["F3"]),
                "F4"      = unname(fibrosis_stage_colors["F4"]))

stage_n <- scores_stage[, .(N = .N, med = round(median(mean_lfc_up), 3)),
                         by = stage_group]
message("Mean log2FC up-DEGs by stage group:"); print(stage_n[order(stage_group)])

pE1 <- ggplot(scores_stage, aes(x = stage_group, y = mean_lfc_up, fill = stage_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.45) +
  geom_violin(alpha = 0.5, color = "white", trim = TRUE, scale = "width") +
  geom_jitter(width = 0.10, height = 0, size = 0.5, alpha = 0.32, color = "grey25") +
  geom_boxplot(width = 0.12, outlier.shape = NA, color = "grey15", fill = "white",
               alpha = 0.5, linewidth = 0.42) +
  geom_text(data = stage_n, aes(x = stage_group, y = -Inf, label = paste0("n=", N)),
            vjust = -0.3, size = 2.6, color = "grey45", inherit.aes = FALSE) +
  scale_fill_manual(values = stage_fill, guide = "none") +
  scale_x_discrete(labels = c("Control" = "Control\n(healthy)",
                               "F0" = "F0\n(disease)", "F1" = "F1",
                               "F2" = "F2", "F3" = "F3", "F4" = "F4")) +
  labs(
    x        = "Group  (healthy controls → fibrosis stage)",
    y        = expression("Mean log"[2]*"FC per patient"),
    title    = "Per-patient mean log₂FC of dream DEGs across fibrosis stages"
  )

ggsave(file.path(OUTDIR, "E1_fibrosis_stage_score.pdf"),
       pE1, width = 7, height = 5, device = cairo_pdf)
message("Saved: E1_fibrosis_stage_score.pdf")

# ── E1 (MEDIAN variant, for comparison vs the mean) ──────────────────────────
# Same sign-oriented per-gene LFC (up as-is, down flipped), but per patient take the
# MEDIAN over DEGs instead of the mean. Robust to outlier genes; compresses dynamic range.
med_score <- function(mat, up, dn) {
  oriented <- rbind(
    if (length(up)) mat[up, , drop = FALSE]  else NULL,
    if (length(dn)) -mat[dn, , drop = FALSE] else NULL)
  apply(oriented, 2, median, na.rm = TRUE)
}
med_dis  <- med_score(lfc_m,    up_g, dn_g)
med_ctrl <- med_score(ctrl_lfc, up_g, dn_g)
scores_dis[,  med_lfc := med_dis[sample_id]]
scores_ctrl[, med_lfc := med_ctrl[sample_id]]
message(sprintf("  Median log2FC disease direction — disease median=%.3f, control median=%.3f",
                median(med_dis), median(med_ctrl)))
message(sprintf("  CONTROL spread — mean-metric   [min,Q1,med,Q3,max]: %s",
                paste(round(quantile(mean_lfc_ctrl, c(0,.25,.5,.75,1)), 4), collapse=" / ")))
message(sprintf("  CONTROL spread — median-metric [min,Q1,med,Q3,max]: %s",
                paste(round(quantile(med_ctrl, c(0,.25,.5,.75,1)), 4), collapse=" / ")))
message(sprintf("  mean of control scores (mean-metric) = %.3e  (algebraically ~0 by construction)",
                mean(mean_lfc_ctrl)))
message(sprintf("  controls landing exactly on 0 (median-metric): %d / %d",
                sum(med_ctrl == 0), length(med_ctrl)))

scores_stage_med <- copy(rbind(scores_dis, scores_ctrl))
scores_stage_med[, stage_group := ifelse(
  group_binary == "Control", "Control",
  ifelse(fibrosis_stage %in% as.character(0:4), paste0("F", fibrosis_stage), NA_character_))]
scores_stage_med <- scores_stage_med[!is.na(stage_group)]
scores_stage_med[, stage_group := factor(stage_group, levels = c("Control", paste0("F", 0:4)))]
stage_n_med <- scores_stage_med[, .(N = .N), by = stage_group]

pE1_med <- ggplot(scores_stage_med, aes(x = stage_group, y = med_lfc, fill = stage_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.45) +
  geom_violin(alpha = 0.5, color = "white", trim = TRUE, scale = "width") +
  geom_jitter(width = 0.10, height = 0, size = 0.5, alpha = 0.32, color = "grey25") +
  geom_boxplot(width = 0.12, outlier.shape = NA, color = "grey15", fill = "white",
               alpha = 0.5, linewidth = 0.42) +
  geom_text(data = stage_n_med, aes(x = stage_group, y = -Inf, label = paste0("n=", N)),
            vjust = -0.3, size = 2.6, color = "grey45", inherit.aes = FALSE) +
  scale_fill_manual(values = stage_fill, guide = "none") +
  scale_x_discrete(labels = c("Control" = "Control\n(healthy)",
                               "F0" = "F0\n(disease)", "F1" = "F1",
                               "F2" = "F2", "F3" = "F3", "F4" = "F4")) +
  labs(
    x        = "Group  (healthy controls → fibrosis stage)",
    y        = expression("Median log"[2]*"FC per patient"),
    title    = "Per-patient median log₂FC of dream DEGs across fibrosis stages"
  )

ggsave(file.path(OUTDIR, "E1_fibrosis_stage_score_median.pdf"),
       pE1_med, width = 7, height = 5, device = cairo_pdf)
message("Saved: E1_fibrosis_stage_score_median.pdf")

# ════════════════════════════════════════════════════════════════════════════
# E2: Waterfall — every sample ranked by signature score, colored by group
# ════════════════════════════════════════════════════════════════════════════
message("\n── E2: Waterfall ──")

scores_wf <- copy(scores_all)[order(sig_score)]
scores_wf[, rank := .I]

# Color bars: control = gray; disease = fibrosis-stage gradient, or teal if no stage
scores_wf[, bar_fill := col_control]
scores_wf[group_binary == "Disease", bar_fill := col_disease]
for (fs in 0:4) {
  fkey <- paste0("F", fs)
  scores_wf[group_binary == "Disease" & fibrosis_stage == as.character(fs),
            bar_fill := unname(fibrosis_stage_colors[fkey])]
}

n_dis       <- sum(scores_wf$group_binary == "Disease")
n_ctrl      <- sum(scores_wf$group_binary == "Control")
mean_ctrl_r <- mean(scores_wf[group_binary == "Control", rank])
mean_dis_r  <- mean(scores_wf[group_binary == "Disease",  rank])
score_max   <- max(scores_wf$sig_score)
score_95    <- quantile(scores_wf$sig_score, 0.98)

pE2 <- ggplot(scores_wf, aes(x = rank, y = sig_score, fill = bar_fill)) +
  geom_col(width = 1, color = NA, alpha = 0.88) +
  geom_hline(yintercept = 0, color = "grey35", linewidth = 0.35) +
  annotate("text", x = mean_ctrl_r, y = score_95 * 0.88,
           label = sprintf("Controls\n(n=%d)", n_ctrl),
           color = "grey45", size = 3, hjust = 0.5, fontface = "bold") +
  annotate("text", x = mean_dis_r, y = score_95 * 0.88,
           label = sprintf("Disease\n(n=%d)", n_dis),
           color = col_disease, size = 3, hjust = 0.5, fontface = "bold") +
  scale_fill_identity() +
  scale_x_continuous(expand = c(0.005, 0)) +
  labs(
    x        = "All samples ranked by signature score (low → high)",
    y        = "Signature score",
    title    = "Waterfall: healthy controls cluster at 0, disease patients score higher",
    subtitle = sprintf("padj < 0.05, |log₂FC| > %.1f", CUT)
  ) +
  theme(axis.text.x  = element_blank(),
        axis.ticks.x = element_blank())

ggsave(file.path(OUTDIR, "E2_waterfall.pdf"),
       pE2, width = 9, height = 4.5, device = cairo_pdf)
message("Saved: E2_waterfall.pdf")

# ════════════════════════════════════════════════════════════════════════════
# E3: ROC curve — individual-level Disease vs Control discrimination
# ════════════════════════════════════════════════════════════════════════════
message("\n── E3: ROC curve ──")

roc_dt   <- scores_all[, .(sig_score, y = as.integer(group_binary == "Disease"))]
roc_dt   <- roc_dt[order(-sig_score)]
n_pos    <- sum(roc_dt$y == 1)
n_neg    <- sum(roc_dt$y == 0)

# Trapezoidal ROC (no external package)
tp <- c(0, cumsum(roc_dt$y == 1)) / n_pos
fp <- c(0, cumsum(roc_dt$y == 0)) / n_neg
auroc <- sum(diff(fp) * (head(tp, -1) + tail(tp, -1)) / 2)
message(sprintf("  AUROC = %.4f  (n_disease=%d, n_control=%d)", auroc, n_pos, n_neg))

# Youden's J optimal point
j_idx   <- which.max(tp[-1] - fp[-1])
opt_tpr <- tp[j_idx + 1]
opt_fpr <- fp[j_idx + 1]
opt_thr <- roc_dt$sig_score[j_idx]
message(sprintf("  Optimal threshold (Youden J): %.3f  TPR=%.3f  FPR=%.3f",
                opt_thr, opt_tpr, opt_fpr))

roc_curve <- data.table(fpr = fp, tpr = tp)

pE3 <- ggplot(roc_curve, aes(x = fpr, y = tpr)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",
              color = "grey70", linewidth = 0.45) +
  geom_path(color = col_disease, linewidth = 1.15) +
  geom_point(aes(x = opt_fpr, y = opt_tpr),
             color = "#C9265E", size = 2.8, shape = 21, fill = "white", stroke = 1.5) +
  annotate("label",
           x = 0.62, y = 0.20,
           label = sprintf("AUROC = %.3f", auroc),
           size = 4.5, fontface = "bold", color = col_disease,
           fill = "white", label.size = 0.3) +
  annotate("text",
           x = opt_fpr + 0.04, y = opt_tpr - 0.05,
           label = sprintf("Youden J\n(score = %.2f)", opt_thr),
           size = 2.7, color = "#C9265E", hjust = 0) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                     labels = paste0(seq(0, 100, 20), "%")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                     labels = paste0(seq(0, 100, 20), "%")) +
  labs(
    x        = "False positive rate  (1 − specificity)",
    y        = "True positive rate  (sensitivity)",
    title    = "Individual-level discrimination: signature score classifies disease vs. healthy",
    subtitle = sprintf("padj < 0.05, |log₂FC| > %.1f", CUT)
  )

ggsave(file.path(OUTDIR, "E3_roc_curve.pdf"),
       pE3, width = 5.5, height = 5.5, device = cairo_pdf)
message("Saved: E3_roc_curve.pdf")

# ════════════════════════════════════════════════════════════════════════════
# E4: Per-cohort violin + boxplot — same story in every canonical cohort
# ════════════════════════════════════════════════════════════════════════════
message("\n── E4: Per-cohort boxplots ──")

scores_cohort <- scores_all[dataset %in% valid_ds & !is.na(dataset)]

# Per-cohort n and Wilcoxon p
cohort_stats <- scores_cohort[, {
  d <- sig_score[group_binary == "Disease"]
  c <- sig_score[group_binary == "Control"]
  wt <- suppressWarnings(wilcox.test(d, c, exact = FALSE))
  .(n_dis   = length(d),
    n_ctrl  = length(c),
    med_dis = round(median(d), 3),
    med_ctrl = round(median(c), 3),
    p_wil   = wt$p.value)
}, by = dataset]
cohort_stats[, p_label := ifelse(p_wil < 0.001, "p < 0.001",
                          ifelse(p_wil < 0.01,  sprintf("p = %.3f", p_wil),
                                                 sprintf("p = %.2f", p_wil)))]
message("Per-cohort summary:"); print(cohort_stats)

# y-position for p-label: just above the 99th percentile per cohort
cohort_ytop <- scores_cohort[, .(y_top = quantile(sig_score, 0.99)), by = dataset]
cohort_stats <- merge(cohort_stats, cohort_ytop, by = "dataset")

# Cohort x-label: add n per group below name
scores_cohort[, ds_label := dataset]

pE4 <- ggplot(scores_cohort,
              aes(x = group_binary, y = sig_score, fill = group_binary)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.35) +
  geom_violin(alpha = 0.60, color = "white", trim = TRUE, scale = "width") +
  geom_boxplot(width = 0.14, outlier.size = 0.35, outlier.alpha = 0.35,
               color = "grey30", fill = "white", linewidth = 0.40) +
  geom_text(data = cohort_stats,
            aes(x = 1.5, y = y_top * 1.02, label = p_label),
            size = 2.6, inherit.aes = FALSE, color = "grey35", hjust = 0.5) +
  facet_wrap(~ dataset, nrow = 1) +
  scale_fill_manual(values = grp_colors, name = NULL) +
  scale_x_discrete(labels = c("Control" = "Ctrl", "Disease" = "Disease")) +
  labs(
    x        = NULL,
    y        = "Signature score",
    title    = "Consistent separation across all 5 canonical cohorts",
    subtitle = sprintf("padj < 0.05, |log₂FC| > %.1f", CUT)
  ) +
  theme(legend.position = "none",
        axis.text.x     = element_text(size = 8),
        strip.text      = element_text(size = 8.5))

ggsave(file.path(OUTDIR, "E4_cohort_boxplots.pdf"),
       pE4, width = 11, height = 5, device = cairo_pdf)
message("Saved: E4_cohort_boxplots.pdf")

message("\n========== All E-panels saved to: ", OUTDIR, " ==========")
