#!/usr/bin/env Rscript
# ============================================================================
# 349_stage_ccc_figures.R
#
# Stage-trajectory CCC figure: heatmap of top stage-progressive LR pairs
# (rows) x 4 disease-stage bins (cols), with right-margin annotations:
#   - continuous pseudotime beta (from Script 346 continuous axis)
#   - documented F-stage trajectory (if available)
#   - bulk concordance flag (from Script 348)
#   - LOO replication rate (from Script 347)
#
# Output: figures/main/fig3_RNAseq/panels/figS_stage_ccc_trajectory.pdf
# (placed in supplementary by default; promote to figs/main if user approves)
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUT_DIR <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory")
FIG_DIR <- file.path(BASE, "figures/main/fig3_RNAseq/panels")
SUPP_DIR<- file.path(BASE, "figures/supplementary/stage_ccc")
dir.create(SUPP_DIR, showWarnings = FALSE, recursive = TRUE)

source(file.path(BASE, "scripts/figures/publication_theme.R"))

# Inputs --------------------------------------------------------------------
COARSE   <- fread(file.path(OUT_DIR, "stage_lr_lmm_coarse.tsv"))
META_EXT <- fread(file.path(OUT_DIR, "donor_metadata_extended.tsv"))
MERGED_TSV <- file.path(OUT_DIR, "all_donor_lr_scores.tsv.gz")

# Optional inputs
CONTINUOUS <- if (file.exists(file.path(OUT_DIR, "stage_lr_lmm_continuous.tsv"))) {
  fread(file.path(OUT_DIR, "stage_lr_lmm_continuous.tsv"))
} else NULL
BULK_CONC  <- if (file.exists(file.path(OUT_DIR, "lr_bulk_concordance.tsv"))) {
  fread(file.path(OUT_DIR, "lr_bulk_concordance.tsv"))
} else NULL
LOO_RATE   <- if (file.exists(file.path(OUT_DIR, "loo_replication_rate_per_lr.tsv"))) {
  fread(file.path(OUT_DIR, "loo_replication_rate_per_lr.tsv"))
} else NULL

# Pick top-N LR pairs by effect size x significance (SH vs Healthy) --------
sh_term <- "disease_stage_coarseSteatohepatitis"
top <- COARSE[term == sh_term]
top[, neglog10p := -log10(pmax(pval, 1e-30))]
top[, rank_score := abs(Estimate) * neglog10p]
top <- top[order(-rank_score)][, head(.SD, 60)]
cat(sprintf("[fig] selected top 60 stage-progressive LR pairs (by |E| x -log10(p))\n"))

# Load merged TSV (produced by Script 345b) --------------------------------
if (!file.exists(MERGED_TSV)) {
  stop("Merged LR scores TSV not found - run Script 345b first")
}
lr_long <- fread(MERGED_TSV)
lr_long[, lr_pair := paste(ligand_complex, receptor_complex, sep = "__")]
lr_long[, ct_pair := paste(source, target, sep = "->")]
lr_long <- merge(lr_long,
                 META_EXT[, .(sample, disease_stage_coarse)],
                 by = "sample", all.x = TRUE)
lr_long[, score := -log10(pmax(magnitude_rank, 1e-4))]

top_keys <- unique(top[, .(ct_pair, lr_pair)])
hm <- merge(lr_long, top_keys, by = c("ct_pair", "lr_pair"))
hm_mean <- hm[, .(score_mean = mean(score, na.rm = TRUE),
                  n_donors   = .N),
              by = .(ct_pair, lr_pair, disease_stage_coarse)]

# Z-score within each LR pair so colour highlights stage differences --------
hm_mean[, z := scale(score_mean)[, 1], by = .(ct_pair, lr_pair)]
hm_mean[, disease_stage_coarse := factor(disease_stage_coarse,
        levels = c("Healthy", "Steatosis", "Steatohepatitis", "Cirrhosis"))]

# Order LR pairs by SH-Healthy effect direction so the heatmap reads as a
# trajectory (up-going at top, down-going at bottom).
top[, pair_id := paste(ct_pair, lr_pair, sep = " | ")]
hm_mean[, pair_id := paste(ct_pair, lr_pair, sep = " | ")]
ord <- top[order(-Estimate)]$pair_id
hm_mean[, pair_id := factor(pair_id, levels = ord)]

# Right-margin track: continuous-axis beta, bulk-concordance flag, LOO rep
annot <- top[, .(pair_id, Estimate)]
setnames(annot, "Estimate", "coarse_SH_beta")
if (!is.null(CONTINUOUS)) {
  cont_sh <- CONTINUOUS[term == "macrophage_pseudotime_mean",
                        .(pair_id = paste(paste(source, target, sep = "->"),
                                          paste(ligand_complex, receptor_complex,
                                                sep = "__"),
                                          sep = " | "),
                          continuous_beta = Estimate,
                          continuous_p = pval)]
  annot <- merge(annot, cont_sh, by = "pair_id", all.x = TRUE)
}
if (!is.null(BULK_CONC)) {
  bc <- unique(BULK_CONC[, .(pair_id = paste(paste(source, target, sep = "->"),
                                             paste(ligand_complex, receptor_complex,
                                                   sep = "__"),
                                             sep = " | "),
                             bulk_both = both_concordant)])
  annot <- merge(annot, bc, by = "pair_id", all.x = TRUE)
}
if (!is.null(LOO_RATE)) {
  lr <- LOO_RATE[, .(pair_id = paste(ct_pair, lr_pair, sep = " | "),
                     replication_rate)]
  annot <- merge(annot, lr, by = "pair_id", all.x = TRUE)
}

# MASLD project palette (publication_theme.R `masld_colors`):
#   diverging: down (blue) -> white -> up (magenta)
DIV_LOW  <- masld_colors$down   # "#1565C0"
DIV_HIGH <- masld_colors$up     # "#C2185B"

# Main heatmap --------------------------------------------------------------
p_hm <- ggplot(hm_mean,
               aes(x = disease_stage_coarse, y = pair_id, fill = z)) +
  geom_tile(color = "white", linewidth = 0.12) +
  scale_fill_gradient2(low = DIV_LOW, mid = "white", high = DIV_HIGH,
                       midpoint = 0,
                       limits = c(-2, 2), oob = scales::squish,
                       name = "z(-log10 mag.rank)") +
  scale_x_discrete(labels = c("Healthy", "Steat.", "SH", "Cirr.")) +
  labs(x = NULL, y = NULL,
       title = "Stage-trajectory CCC: top LR pairs") +
  theme_masld(base_size = 7) +
  theme(axis.text.y      = element_text(size = 5.5),
        axis.text.x      = element_text(angle = 30, hjust = 1, size = 7),
        axis.ticks       = element_blank(),
        axis.line        = element_blank(),
        legend.position  = "top",
        legend.key.width = unit(0.6, "cm"),
        legend.key.height= unit(0.25, "cm"),
        legend.title     = element_text(size = 6, face = "bold"),
        legend.text      = element_text(size = 5.5),
        legend.margin    = margin(0, 0, 0, 0),
        plot.title       = element_text(face = "bold", size = 9,
                                        hjust = 0, margin = margin(b = 2)),
        plot.title.position = "plot",
        plot.margin      = margin(3, 3, 3, 3))

# Right-margin annotation bars ---------------------------------------------
annot[, pair_id := factor(pair_id, levels = ord)]
annot_long <- melt(annot, id.vars = "pair_id",
                   measure.vars = intersect(c("continuous_beta",
                                               "bulk_both",
                                               "replication_rate"),
                                            names(annot)),
                   variable.name = "track", value.name = "value")
annot_long[, value := as.numeric(value)]
annot_long[, track := factor(track,
  levels = c("continuous_beta", "bulk_both", "replication_rate"),
  labels = c("PT β", "Bulk", "LOO"))]

p_annot <- ggplot(annot_long,
                  aes(x = track, y = pair_id, fill = value)) +
  geom_tile(color = "white", linewidth = 0.12) +
  scale_fill_gradient2(low = DIV_LOW, mid = "white", high = DIV_HIGH,
                       midpoint = 0, na.value = "grey90",
                       name = "value") +
  labs(x = NULL, y = NULL, title = "Tracks") +
  theme_masld(base_size = 6) +
  theme(axis.text.y      = element_blank(),
        axis.text.x      = element_text(angle = 30, hjust = 1, size = 6),
        axis.ticks       = element_blank(),
        axis.line        = element_blank(),
        legend.position  = "top",
        legend.key.width = unit(0.4, "cm"),
        legend.key.height= unit(0.2, "cm"),
        legend.title     = element_text(size = 5.5, face = "bold"),
        legend.text      = element_text(size = 5),
        legend.margin    = margin(0, 0, 0, 0),
        plot.title       = element_text(face = "bold", size = 7,
                                        hjust = 0, margin = margin(b = 2)),
        plot.title.position = "plot",
        plot.margin      = margin(3, 3, 3, 3))

panel <- (p_hm | p_annot) + plot_layout(widths = c(4, 1))

out_pdf <- file.path(SUPP_DIR, "figS_stage_ccc_trajectory.pdf")
save_fig(panel, out_pdf, width = 8, height = 8)
cat(sprintf("[output] %s\n", out_pdf))

# Companion: per-cell-type-pair count of emergent vs declining LR pairs ----
ct_summary <- top[, .(n_emergent = sum(Estimate > 0 & pval < 0.05),
                      n_declining = sum(Estimate < 0 & pval < 0.05)),
                  by = ct_pair]
ct_summary[, total := n_emergent + n_declining]
ct_summary <- ct_summary[order(-total)][, head(.SD, 20)]
ct_long <- melt(ct_summary, id.vars = c("ct_pair", "total"),
                measure.vars = c("n_emergent", "n_declining"),
                variable.name = "direction", value.name = "n")

p_ct <- ggplot(ct_long, aes(x = reorder(ct_pair, total), y = n, fill = direction)) +
  geom_col(width = 0.72) +
  coord_flip() +
  scale_fill_manual(values = c(n_emergent  = masld_colors$up,
                               n_declining = masld_colors$down),
                    labels = c("Emergent (up in SH)",
                               "Declining (down in SH)"),
                    name = NULL) +
  labs(x = NULL, y = "# significant LR pairs (top 60)",
       title = "Top 20 cell-type pairs by stage-progressive LR count") +
  theme_masld(base_size = 7) +
  theme(legend.position    = "top",
        legend.key.size    = unit(0.3, "cm"),
        legend.text        = element_text(size = 6),
        legend.margin      = margin(0, 0, 0, 0),
        legend.box.spacing = unit(0.1, "cm"),
        axis.text.y        = element_text(size = 6),
        axis.text.x        = element_text(size = 6),
        axis.title.x       = element_text(size = 7),
        plot.title         = element_text(face = "bold", size = 9,
                                          hjust = 0, margin = margin(b = 3)),
        plot.title.position= "plot",
        plot.margin        = margin(4, 5, 4, 4))

out_pdf2 <- file.path(SUPP_DIR, "figS_stage_ccc_ct_summary.pdf")
save_fig(p_ct, out_pdf2, width = 6, height = 5)
cat(sprintf("[output] %s\n", out_pdf2))

cat("\n[done] figures saved to figures/supplementary/stage_ccc/\n")
