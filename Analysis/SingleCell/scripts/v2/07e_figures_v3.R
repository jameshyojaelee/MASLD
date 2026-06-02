#!/usr/bin/env Rscript
# ============================================================================
# 07e_figures_v3.R
#
# Hardened replacement for figS_stage_ccc_trajectory_v2_atlas.pdf.
# Reads stage_lr_headline_v3.tsv (from 07d_chain_v3_final.R) and renders:
#   (i)   heatmap of headline LR pairs (rows) x stage bins (cols), z-scored.
#   (ii)  right-margin strip: documented-F-stage Estimate +/- SE.
#   (iii) right-margin strip: bootstrap 95% CI (Q025 .. Q975).
#   (iv)  optional side panel: three-way method concordance (n_methods_in_topN).
#   (v)   row annotation: endocrine_suspect flag.
#
# Output: figures/supplementary/stage_ccc/figS_stage_ccc_hardened.pdf
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
V3_DIR <- file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2_phase05/ccc/stage_trajectory_v3")
V2_DIR <- file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2_phase05/ccc/stage_trajectory_v2")
FIG_DIR <- file.path(BASE, "figures/supplementary/stage_ccc")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
MERGED_TSV <- file.path(V2_DIR, "all_donor_lr_scores_v2.tsv.gz")
HEADLINE_TSV <- file.path(V3_DIR, "stage_lr_headline_v3.tsv")

if (!file.exists(HEADLINE_TSV)) stop("missing headline list at ", HEADLINE_TSV)

hl <- fread(HEADLINE_TSV)
cat(sprintf("[input] %d total LR pairs in headline table\n", nrow(hl)))
cat(sprintf("[input] %d pass all gates (headline=TRUE)\n",
            sum(hl$headline, na.rm = TRUE)))

# Pick top-60 by rank_score from the FULL table (not just headline-pass)
# so that the heatmap always has 60 rows. Mark which rows pass.
hl[, neg_rs := -rank_score]
setorder(hl, neg_rs)
top60 <- hl[seq_len(min(60, nrow(hl)))]
top60[, pair_id := paste(ct_pair, lr_pair, sep = " | ")]

# Load donor-LR scores to compute per-stage means
lr_long <- fread(MERGED_TSV)
lr_long[, lr_pair := paste(ligand_complex, receptor_complex, sep = "__")]
lr_long[, ct_pair := paste(source, target, sep = "->")]
lr_long[, score := -log10(pmax(magnitude_rank, 1e-4))]
lr_long[, pair_id := paste(ct_pair, lr_pair, sep = " | ")]

# Add donor stage from v2 metadata
META_V2 <- file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2_phase05/mcp/inputs/donor_metadata_v2.tsv")
meta <- fread(META_V2)
stage_levels <- c("Healthy", "Steatosis", "Steatohepatitis", "Cirrhosis")
meta_use <- meta[, .(sample, disease_stage_coarse)]
meta_use <- meta_use[!is.na(disease_stage_coarse) & disease_stage_coarse != ""]
meta_use[, disease_stage_coarse := factor(disease_stage_coarse, levels = stage_levels)]

lr_long <- merge(lr_long, meta_use, by = "sample", all.x = FALSE)
hm <- lr_long[pair_id %in% top60$pair_id]
hm_mean <- hm[, .(score_mean = mean(score, na.rm = TRUE),
                  n_donors = .N),
              by = .(pair_id, disease_stage_coarse)]
hm_mean[, z := scale(score_mean)[, 1], by = pair_id]
hm_mean[, pair_id := factor(pair_id, levels = top60$pair_id)]

# --- Shared label formatter: abbreviated cell types + "→" separators -----
fmt_pair <- function(x) {
  x <- gsub("->",               "→",     x, fixed = TRUE)
  x <- gsub("__",               "→",     x, fixed = TRUE)
  x <- gsub("Endothelial cells","Endo",  x)
  x <- gsub("Hepatocytes",      "Hep",   x)
  x <- gsub("Macrophages",      "Mac",   x)
  x <- gsub("Fibroblasts",      "Fib",   x)
  x <- gsub("Cholangiocytes",   "Chol",  x)
  x <- gsub("T cells",          "Tcell", x)
  x <- gsub("B cells",          "Bcell", x)
  x
}

# --- Panel A: stage z-heatmap -------------------------------------------
hm_mean[, pair_label := fmt_pair(as.character(pair_id))]
hm_mean[, pair_label := factor(pair_label,
  levels = fmt_pair(levels(hm_mean$pair_id)))]

p_heat <- ggplot(hm_mean, aes(x = disease_stage_coarse, y = pair_label, fill = z)) +
  geom_tile(color = "white", linewidth = 0.15) +
  scale_fill_gradient2(low = masld_colors$down, mid = "#FFFFFF", high = masld_colors$up,
                       midpoint = 0, name = "z-score",
                       limits = c(-1.5, 1.5), oob = scales::squish) +
  labs(x = NULL, y = NULL,
       title = "Top 60 stage-progressive LR pairs",
       subtitle = "Mean –log10(magnitude rank) z-scored per pair") +
  theme_masld(base_size = 7) + theme_pub() +
  theme(axis.text.y   = element_text(size = 4.5, lineheight = 0.9),
        axis.text.x   = element_text(angle = 30, hjust = 1),
        legend.position = "bottom",
        legend.key.width = unit(0.4, "cm"),
        legend.key.height = unit(0.12, "cm"))

# --- Panel B: documented-F-stage estimate +/- SE ------------------------
if ("doc_Estimate" %in% names(top60) && any(!is.na(top60$doc_Estimate))) {
  doc_dt <- top60[, .(pair_id = factor(pair_id, levels = top60$pair_id),
                      doc_Estimate, doc_pval)]
  doc_dt[, pair_label := factor(fmt_pair(as.character(pair_id)),
                                levels = fmt_pair(levels(pair_id)))]
  p_doc <- ggplot(doc_dt, aes(x = doc_Estimate, y = pair_label)) +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey70", linetype = "dashed") +
    geom_point(aes(color = doc_Estimate > 0), size = 0.9, stroke = 0) +
    scale_color_manual(values = c(`TRUE` = masld_colors$up, `FALSE` = masld_colors$down),
                       guide = "none", na.value = "grey60") +
    scale_x_continuous(name = "F-stage β", n.breaks = 3) +
    theme_masld(base_size = 7) + theme_pub() +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank()) +
    labs(title = "Documented\nF-stage")
} else {
  p_doc <- ggplot() + ggtitle("Documented F-stage\n(no data)") +
    theme_void(base_size = 7)
}

# --- Panel C: bootstrap CI ----------------------------------------------
if ("boot_median" %in% names(top60) && any(!is.na(top60$boot_median))) {
  boot_dt <- top60[, .(pair_id = factor(pair_id, levels = top60$pair_id),
                       boot_median, boot_q025, boot_q975, ci_excludes_zero)]
  boot_dt[, pair_label := factor(fmt_pair(as.character(pair_id)),
                                 levels = fmt_pair(levels(pair_id)))]
  p_boot <- ggplot(boot_dt, aes(y = pair_label)) +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey70", linetype = "dashed") +
    geom_errorbarh(aes(xmin = boot_q025, xmax = boot_q975,
                       color = ci_excludes_zero),
                   height = 0, linewidth = 0.4) +
    geom_point(aes(x = boot_median, color = ci_excludes_zero), size = 0.8) +
    scale_color_manual(values = c("TRUE"  = masld_colors$nash,
                                  "FALSE" = masld_colors$ns),
                       name = "CI excl. 0",
                       labels = c("TRUE" = "Yes", "FALSE" = "No"),
                       na.value = "#BDBDBD") +
    scale_x_continuous(name = "Bootstrap 95% CI", n.breaks = 3) +
    theme_masld(base_size = 7) + theme_pub() +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          legend.position = c(0.78, 0.12),
          legend.background = element_blank()) +
    labs(title = "Bootstrap CI\n(SH vs Healthy)")
} else {
  p_boot <- ggplot() + ggtitle("Bootstrap CI\n(no data)") +
    theme_void(base_size = 7)
}

# --- Panel D: gate / endocrine annotation ------------------------------
ann_cols <- c("gate_significance", "gate_bootstrap", "gate_permutation",
              "gate_leverage", "gate_threeway", "headline",
              "endocrine_suspect")
ann_cols <- intersect(ann_cols, names(top60))
ann_long <- melt(top60[, c("pair_id", ann_cols), with = FALSE],
                 id.vars = "pair_id",
                 measure.vars = ann_cols,
                 variable.name = "gate", value.name = "pass")
ann_long[, pair_id := factor(pair_id, levels = top60$pair_id)]
ann_long[, gate := factor(gate, levels = ann_cols)]
ann_long[, pass_label := ifelse(is.na(pass), "NA",
                                ifelse(pass, "Pass", "Fail"))]

ann_long[, pair_label := factor(fmt_pair(as.character(pair_id)),
                               levels = fmt_pair(levels(ann_long$pair_id)))]
gate_labels <- c(gate_significance = "Signif.",
                 gate_bootstrap    = "Bootstrap",
                 gate_permutation  = "Permut.",
                 gate_leverage     = "Leverage",
                 gate_threeway     = "3-way",
                 headline          = "Headline",
                 endocrine_suspect = "Endocrine?")
ann_long[, gate_label := factor(gate_labels[as.character(gate)],
                                levels = gate_labels)]

p_ann <- ggplot(ann_long, aes(x = gate_label, y = pair_label, fill = pass_label)) +
  geom_tile(color = "white", linewidth = 0.15) +
  scale_fill_manual(values = c(Pass = masld_colors$nash,
                               Fail = "#E0E0E0",
                               `NA` = "#F5F5F5"),
                    name = NULL) +
  theme_masld(base_size = 7) + theme_pub() +
  theme(axis.text.y  = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x  = element_text(size = 5, angle = 40, hjust = 1, vjust = 1),
        axis.line    = element_blank(),
        axis.ticks   = element_blank()) +
  labs(x = NULL, y = NULL, title = "Quality gates")

# --- Compose ------------------------------------------------------------
pdf_path <- file.path(FIG_DIR, "figS_stage_ccc_hardened.pdf")
combined <- (p_heat | p_doc | p_boot | p_ann) +
  plot_layout(widths = c(2.2, 0.8, 1.2, 1.4)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(size = 8, face = "bold",
                                          margin = margin(0, 2, 0, 0)))
  )
ggsave(pdf_path, combined, width = 11, height = 11, device = cairo_pdf)
cat(sprintf("[output] %s\n", pdf_path))
cat("[done]\n")
