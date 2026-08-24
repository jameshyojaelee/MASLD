#!/usr/bin/env Rscript
##############################################################################
# fig5_panel_calibration.R -- Fig 5 panel 5d (was 5e; renamed 2026-07-08)
#
# KEY MESSAGE: External MASLD benchmark panels are recovered by different
# evidence channels; traditional ROC curves show both benchmark recovery and
# the limits of treating convergence as a single-channel replacement.
#
# This panel is an external-panel recovery check, not a calibration curve and
# not evidence that integration maximizes AUROC. It deliberately uses one
# evaluation universe and one label vector per panel for every plotted score.
#
# HONEST FRAMING (robustness, NOT superiority): each single channel wins its
# own-modality benchmark (Expression -> expression signatures; Genetics ->
# genetic loci) and degrades off-modality (Expression falls to ~0.50 on genetic
# loci). Full-stack convergence is never panel-best on the genetic/clinical
# benchmarks, but has the highest mean and highest worst-case AUROC across the
# five modalities -- it generalizes where no single channel does. Per-panel
# AUROCs are annotated in-facet; head-to-head DeLong tests (convergence vs each
# channel) are written to fig6_calibration_delong_pairwise.csv.
##############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

set.seed(42)
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANDIR <- file.path(FIG5_DIR, "panels")
dir.create(PANDIR, recursive = TRUE, showWarnings = FALSE)

ME_DIR <- file.path(BASE, "RNA-seq/results/multi_evidence")
PANELS_DIR <- file.path(BASE, "data/published_gene_panels")

# Known legacy/public-resource aliases that otherwise create false missing genes.
ALIAS_MAP <- c(
  PPAPDC1A = "PLPP4",
  ACC1     = "ACACA",
  FXR      = "NR1H4",
  ASK1     = "MAP3K5",
  GAL3     = "LGALS3",
  VAP1     = "AOC3",
  SGLT2    = "SLC5A2",
  MARC1    = "MTARC1"
)

canonicalize_genes <- function(x) {
  x <- unique(trimws(x))
  x <- x[!is.na(x) & nzchar(x)]
  mapped <- ifelse(x %chin% names(ALIAS_MAP), unname(ALIAS_MAP[x]), x)
  unique(mapped)
}

read_panel_genes <- function(path, filter_expr = NULL) {
  raw <- readLines(path)
  raw <- raw[!grepl("^#", raw)]
  d <- fread(text = paste(raw, collapse = "\n"))
  if (!is.null(filter_expr)) d <- d[eval(filter_expr)]
  canonicalize_genes(d$gene_symbol)
}

read_panel_gene_column <- function(path, gene_col, filter_expr = NULL,
                                   split_slash = FALSE) {
  raw <- readLines(path)
  raw <- raw[!grepl("^#", raw)]
  d <- fread(text = paste(raw, collapse = "\n"))
  if (!is.null(filter_expr)) d <- d[eval(filter_expr)]
  genes <- d[[gene_col]]
  if (split_slash) genes <- unlist(strsplit(genes, "/", fixed = TRUE))
  canonicalize_genes(genes)
}

auroc_from_scores <- function(scores, labels) {
  labels <- as.logical(labels)
  r <- rank(scores, ties.method = "average")
  n1 <- sum(labels)
  n0 <- sum(!labels)
  (sum(r[labels]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

hanley_ci <- function(a, n1, n0) {
  if (!is.finite(a) || n1 < 2 || n0 < 2) return(c(lo = NA_real_, hi = NA_real_))
  q1 <- a / (2 - a)
  q2 <- 2 * a^2 / (1 + a)
  se <- sqrt((a * (1 - a) + (n1 - 1) * (q1 - a^2) +
                (n0 - 1) * (q2 - a^2)) / (n1 * n0))
  c(lo = max(0, a - 1.96 * se), hi = min(1, a + 1.96 * se))
}

pr_auc_manual <- function(scores, labels) {
  labels <- as.logical(labels)
  ok <- is.finite(scores) & !is.na(labels)
  scores <- scores[ok]
  labels <- labels[ok]
  if (sum(labels) < 2 || sum(!labels) < 2) return(NA_real_)

  ord <- order(scores, decreasing = TRUE)
  labels <- labels[ord]
  tp <- cumsum(labels)
  fp <- cumsum(!labels)
  precision <- tp / (tp + fp)
  recall <- tp / sum(labels)
  recall <- c(0, recall)
  precision <- c(1, precision)
  sum(diff(recall) * (precision[-1] + precision[-length(precision)]) / 2)
}

roc_curve <- function(scores, labels, max_points = 1000L) {
  labels <- as.logical(labels)
  ok <- is.finite(scores) & !is.na(labels)
  scores <- scores[ok]
  labels <- labels[ok]
  n_pos <- sum(labels)
  n_neg <- sum(!labels)
  if (n_pos < 2 || n_neg < 2) {
    return(data.table(fpr = numeric(), tpr = numeric()))
  }

  d <- data.table(score = scores, label = labels)
  d <- d[, .(tp = sum(label), fp = sum(!label)), by = score]
  setorder(d, -score)
  d[, `:=`(tpr = cumsum(tp) / n_pos, fpr = cumsum(fp) / n_neg)]
  curve <- rbind(
    data.table(fpr = 0, tpr = 0),
    d[!(fpr == 1 & tpr == 1), .(fpr, tpr)],
    data.table(fpr = 1, tpr = 1)
  )
  if (nrow(curve) > max_points) {
    keep <- sort(unique(c(
      1L,
      nrow(curve),
      as.integer(round(seq(1, nrow(curve), length.out = max_points)))
    )))
    curve <- curve[keep]
  }

  curve
}

score_one <- function(dt, panel_name, panel_label, panel_type, panel_genes,
                      score_id, score_label, score_family, scores) {
  valid <- is.finite(scores)
  labels <- dt$human_symbol %chin% panel_genes
  n_pos <- sum(labels[valid])
  n_neg <- sum(!labels[valid])
  if (n_pos < 2 || n_neg < 2) {
    au <- NA_real_
    ci <- c(lo = NA_real_, hi = NA_real_)
    p <- NA_real_
    pr <- NA_real_
  } else {
    au <- auroc_from_scores(scores[valid], labels[valid])
    ci <- hanley_ci(au, n_pos, n_neg)
    p <- suppressWarnings(
      wilcox.test(scores[valid][labels[valid]], scores[valid][!labels[valid]],
                  alternative = "greater", exact = FALSE)$p.value
    )
    pr <- pr_auc_manual(scores[valid], labels[valid])
  }
  missing <- setdiff(panel_genes, dt$human_symbol[valid])
  data.table(
    panel = panel_name,
    label = panel_label,
    panel_type = panel_type,
    n_panel_canonical = length(panel_genes),
    n_pos = n_pos,
    n_neg = n_neg,
    n_missing = length(missing),
    missing_genes = paste(missing, collapse = ";"),
    score_id = score_id,
    score_label = score_label,
    score_family = score_family,
    auroc = au,
    ci_lo = unname(ci["lo"]),
    ci_hi = unname(ci["hi"]),
    ci_method = "Hanley-McNeil",
    pr_auc = pr,
    prevalence = n_pos / (n_pos + n_neg),
    p = p,
    star = fifelse(is.na(p), "NA",
      fifelse(p < 1e-3, "***",
        fifelse(p < 1e-2, "**", fifelse(p < 0.05, "*", "ns"))))
  )
}

# Load score universe.
atlas <- fread(file.path(ME_DIR, "multi_evidence_atlas.csv"),
  select = c("human_symbol", "bulk_logFC", "bulk_padj", "bulk_treat_fdr",
             "coloc_susie_best_pp4", "coloc_abf_best_pp4"))
ce <- fread(file.path(ME_DIR, "convergence_evidence.csv"),
  select = c("human_symbol", "excluded_from_ranking", "convergence_score",
             "log_BF_S1", "log_BF_S2_coloc", "log_BF_S3", "log_BF_S4",
             "log_BF_S5", "log_BF_S6", "log_BF_S7", "log_BF_S8"))

dt <- merge(ce, atlas, by = "human_symbol", all.x = TRUE)
dt <- dt[excluded_from_ranking == FALSE & is.finite(convergence_score)]

bf_cols <- c("log_BF_S2_coloc", "log_BF_S3", "log_BF_S4", "log_BF_S5",
             "log_BF_S6", "log_BF_S7", "log_BF_S8")
dt[, score_convergence_full := convergence_score]
for (cc in bf_cols) dt[is.na(get(cc)) | get(cc) < 0, (cc) := 0]
dt[, score_convergence_lomo_s1 := rowSums(.SD), .SDcols = bf_cols]

dt[, score_top_logFC := fifelse(is_canonical_deg(dt),
                                abs(bulk_logFC), 0)]
dt[is.na(score_top_logFC), score_top_logFC := 0]
dt[, score_top_pp4 := pmax(coloc_susie_best_pp4, coloc_abf_best_pp4, na.rm = TRUE)]
dt[is.na(score_top_pp4) | is.infinite(score_top_pp4), score_top_pp4 := 0]

protein_marker_genes <- canonicalize_genes(c(
  read_panel_genes(
    file.path(PANELS_DIR, "masld_nontranscriptomic_markers.tsv"),
    quote(marker_class %chin% c("plasma_protein", "liver_protein") &
            !is.na(gene_symbol) & gene_symbol != "NA")
  ),
  # Named markers from Govaere/Anstee 2023 Nat Metab proteo-transcriptomic
  # analysis: eight shared NAS/fibrosis markers, named advanced-fibrosis/NAS
  # markers, and cell-origin examples from the 31-marker signature.
  c("THBS2", "APOF", "ADAMTSL2", "CFHR4", "TREM2", "AKR1B10",
    "SULT2A1", "PTGR1", "GDF15", "IGFBP7", "SHBG", "ADSSL1",
    "ENO3", "CXCL8")
))

panels <- list(
  Govaere = list(
    genes = read_panel_genes(file.path(PANELS_DIR, "govaere_2020_panel.tsv")),
    label = "Govaere signature",
    type = "expression"
  ),
  Moylan = list(
    genes = read_panel_genes(file.path(PANELS_DIR, "moylan_2014_panel.tsv")),
    label = "Moylan fibrosis",
    type = "expression"
  ),
  Biomarkers = list(
    genes = protein_marker_genes,
    label = "Protein markers",
    type = "translational"
  ),
  Genetics = list(
    genes = read_panel_gene_column(
      file.path(PANELS_DIR, "mancina_2025_panel.tsv"),
      "gene",
      quote(origin == "germline" &
              stage %chin% c("Steatosis", "Steatohepatitis", "Cirrhosis")),
      split_slash = TRUE
    ),
    label = "MASLD loci",
    type = "genetic"
  ),
  NIDDK = list(
    genes = read_panel_genes(file.path(PANELS_DIR, "niddk_pipeline_2024.tsv")),
    label = "Clinical-trial targets",
    type = "clinical"
  )
)

metric_rows <- list()
for (pn in names(panels)) {
  pinfo <- panels[[pn]]
  metric_rows[[paste(pn, "full", sep = "|")]] <- score_one(
    dt, pn, pinfo$label, pinfo$type, pinfo$genes,
    "score_convergence_full", "Convergence full-stack", "convergence",
    dt$score_convergence_full
  )
  metric_rows[[paste(pn, "lomo_s1", sep = "|")]] <- score_one(
    dt, pn, pinfo$label, pinfo$type, pinfo$genes,
    "score_convergence_lomo_s1", "Convergence LOMO-S1", "convergence",
    dt$score_convergence_lomo_s1
  )
  metric_rows[[paste(pn, "top_logFC", sep = "|")]] <- score_one(
    dt, pn, pinfo$label, pinfo$type, pinfo$genes,
    "score_top_logFC", "top |logFC|", "single_channel", dt$score_top_logFC
  )
  metric_rows[[paste(pn, "top_pp4", sep = "|")]] <- score_one(
    dt, pn, pinfo$label, pinfo$type, pinfo$genes,
    "score_top_pp4", "top PP4", "single_channel", dt$score_top_pp4
  )
}
metrics <- rbindlist(metric_rows, use.names = TRUE)
metrics[, auroc_label := fifelse(star == "ns",
  sprintf("%.2f ns", auroc),
  sprintf("%.2f", auroc))]
fwrite(metrics, file.path(PANDIR, "fig6_calibration_source.csv"))

# Full channel-name lookup (used by the DeLong table); the plotted panel shows
# only the three head-to-head scorers -- the No bulk-DE ablation is retained in
# fig6_calibration_delong_pairwise.csv as a circularity control, not on the main figure.
channel_names <- c(
  score_convergence_full = "Convergence",
  score_convergence_lomo_s1 = "No bulk-DE",
  score_top_logFC = "Expression",
  score_top_pp4 = "Genetics"
)
display_scores <- channel_names[c("score_convergence_full",
                                  "score_top_logFC", "score_top_pp4")]
plot_dt <- metrics[score_id %chin% names(display_scores)]
plot_dt[, channel := unname(display_scores[score_id])]
plot_dt[, channel := factor(channel, levels = unname(display_scores))]
plot_dt[, panel_facet := sprintf("%s\n(%s, n=%d)", label, panel_type, n_pos)]
panel_order <- c("Govaere", "Moylan", "Biomarkers", "Genetics", "NIDDK")
facet_levels <- plot_dt[match(panel_order, panel),
  sprintf("%s\n(%s, n=%d)", label, panel_type, n_pos)]
plot_dt[, panel_facet := factor(panel_facet, levels = facet_levels)]
plot_dt[, is_best_in_panel := auroc == max(auroc, na.rm = TRUE), by = panel]
plot_dt[, best_channel := as.character(channel[which.max(auroc)]), by = panel]
plot_dt[, line_label := sprintf("%s %.2f%s", channel, auroc,
                                fifelse(star == "ns", " ns", ""))]

roc_rows <- list()
score_lookup <- list(
  score_convergence_full = dt$score_convergence_full,
  score_convergence_lomo_s1 = dt$score_convergence_lomo_s1,
  score_top_logFC = dt$score_top_logFC,
  score_top_pp4 = dt$score_top_pp4
)
for (i in seq_len(nrow(plot_dt))) {
  row <- plot_dt[i]
  genes <- panels[[row$panel]]$genes
  labels <- dt$human_symbol %chin% genes
  curve <- roc_curve(score_lookup[[row$score_id]], labels)
  curve[, `:=`(
    panel = row$panel,
    panel_facet = row$panel_facet,
    score_id = row$score_id,
    channel = row$channel,
    line_label = row$line_label,
    auroc = row$auroc,
    star = row$star,
    is_best_in_panel = row$is_best_in_panel
  )]
  roc_rows[[paste(row$panel, row$score_id, sep = "|")]] <- curve
}
roc_dt <- rbindlist(roc_rows, use.names = TRUE)
roc_dt[, panel_facet := factor(panel_facet, levels = facet_levels)]
roc_dt[, channel := factor(channel, levels = unname(display_scores))]
roc_dt[, line_label := factor(line_label,
  levels = plot_dt[order(channel), unique(line_label)])]
fwrite(roc_dt, file.path(PANDIR, "fig6_calibration_baselines.csv"))

## Pairwise DeLong test (correlated ROC, same observations + label vector per
## panel): does full-stack convergence differ from each alternative channel?
## Stars in the panel are vs AUROC=0.5 only; this supplies the head-to-head test.
delong_rows <- list()
for (pn in names(panels)) {
  pinfo <- panels[[pn]]
  labels <- as.integer(dt$human_symbol %chin% pinfo$genes)
  if (sum(labels) < 2 || sum(labels == 0) < 2) next
  roc_conv <- pROC::roc(labels, dt$score_convergence_full,
                        direction = "<", levels = c(0, 1), quiet = TRUE)
  for (sid in c("score_convergence_lomo_s1", "score_top_logFC", "score_top_pp4")) {
    roc_alt <- pROC::roc(labels, dt[[sid]],
                         direction = "<", levels = c(0, 1), quiet = TRUE)
    tt <- pROC::roc.test(roc_conv, roc_alt, method = "delong")
    delong_rows[[paste(pn, sid)]] <- data.table(
      panel = pn, label = pinfo$label, panel_type = pinfo$type,
      comparison = sprintf("Convergence vs %s", unname(channel_names[sid])),
      auroc_convergence = as.numeric(pROC::auc(roc_conv)),
      auroc_alt = as.numeric(pROC::auc(roc_alt)),
      delta_auroc = as.numeric(pROC::auc(roc_conv)) - as.numeric(pROC::auc(roc_alt)),
      delong_p = tt$p.value)
  }
}
delong <- rbindlist(delong_rows)
delong[, delong_fdr := p.adjust(delong_p, "BH")]
fwrite(delong, file.path(PANDIR, "fig6_calibration_delong_pairwise.csv"))

## Per-facet AUROC block (channels black text; legend maps colour/linetype).
auc_annot <- plot_dt[order(panel, channel),
  .(lbl = paste0("AUROC\n",
      paste(sprintf("%s %.2f", channel, auroc), collapse = "\n"))),
  by = .(panel_facet)]

# Visual hierarchy: Convergence is the primary result (saturated colour, thick,
# drawn on top); Expression and Genetics are secondary reference channels shown
# in muted greys (distinguished by lightness + linetype), never competing with
# the Convergence line for attention.
line_cols <- c(
  "Convergence" = masld_colors$mash,
  "Expression" = "#7FA3D1",
  "Genetics" = "#C9A24B"
)
line_types <- c(
  "Convergence" = "solid",
  "Expression" = "42",
  "Genetics" = "longdash"
)
line_widths <- c(
  "Convergence" = 0.85,
  "Expression" = 0.42,
  "Genetics" = 0.42
)

# Draw the two reference channels first, Convergence last (on top).
roc_dt[, .draw := fifelse(channel == "Convergence", 2L, 1L)]
setorder(roc_dt, .draw, panel_facet, channel, fpr)

p5e <- ggplot(roc_dt, aes(x = fpr, y = tpr, color = channel,
                          linetype = channel, linewidth = channel,
                          group = channel)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted",
              linewidth = 0.28, color = "grey60") +
  geom_step(alpha = 0.98) +
  geom_label(data = auc_annot, inherit.aes = FALSE,
             aes(x = 0.985, y = 0.015, label = lbl),
             hjust = 1, vjust = 0, size = 1.75, lineheight = 0.92,
             label.size = 0, label.padding = unit(0.35, "mm"),
             fill = alpha("white", 0.55), color = "black") +
  facet_wrap(~ panel_facet, nrow = 1) +
  scale_color_manual(values = line_cols, name = NULL) +
  scale_linetype_manual(values = line_types, name = NULL) +
  scale_linewidth_manual(values = line_widths, guide = "none") +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1),
                     expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1),
                     expand = expansion(mult = c(0.01, 0.01))) +
  labs(x = "False positive rate", y = "True positive rate") +
  coord_equal() +
  theme_masld(base_size = 6) +
  theme(strip.text = element_text(size = 6, face = "plain", lineheight = 0.8,
                                  margin = margin(b = 1)),
        legend.position = "bottom",
        legend.key.width = unit(0.40, "cm"),
        legend.key.height = unit(0.10, "cm"),
        legend.margin = margin(t = -2),
        legend.box.spacing = unit(0.05, "cm"),
        panel.spacing.x = unit(0.12, "cm"),
        plot.margin = margin(2, 3, 1, 1))

ggsave(file.path(PANDIR, "fig6_calibration.pdf"), p5e,
       width = 6.63, height = 1.88, device = cairo_pdf)

cat("Saved panels/fig6_calibration.pdf + source/baseline CSVs\n\n")
cat("Panel coverage after alias canonicalization:\n")
print(unique(metrics[, .(panel, n_panel_canonical, n_pos, n_missing, missing_genes)]))
cat("\nPlotted ROC summaries:\n")
print(plot_dt[order(panel, channel),
  .(panel, channel, score_label, n_pos, auroc = round(auroc, 3),
    ci = sprintf("[%.2f, %.2f]", ci_lo, ci_hi),
    pr_auc = round(pr_auc, 4), p = signif(p, 2), star, is_best_in_panel)])
cat("\nInterpretation: ROC curves show external benchmark recovery by score channel.\n",
    "Expression-derived benchmarks should favor the expression channel, while genetic loci and\n",
    "clinical targets test whether the genetics/convergence arms recover translational biology.\n")

conv_mean <- plot_dt[, .(mean_auroc = mean(auroc), min_auroc = min(auroc)), by = channel]
cat("\nMean / worst-case AUROC across the five benchmarks (robustness view):\n")
print(conv_mean[order(-mean_auroc)])
cat("\nPairwise DeLong (Convergence full-stack vs each channel):\n")
print(delong[, .(panel, comparison, auroc_convergence = round(auroc_convergence, 3),
  auroc_alt = round(auroc_alt, 3), delta = round(delta_auroc, 3),
  delong_p = signif(delong_p, 2), delong_fdr = signif(delong_fdr, 2))])
