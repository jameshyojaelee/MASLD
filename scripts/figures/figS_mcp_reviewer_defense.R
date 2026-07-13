#!/usr/bin/env Rscript
# figS_mcp_reviewer_defense.R
# Generate reviewer-defense supplementary panels addressing adversarial review:
#   A: cophenetic stability across k (k=10/13/16/20)
#   B: entropy distribution with natural break highlight
#   C: dataset-variance vs stage-effect per program (batch artifact diagnostic)
#   D: changepoint null distribution (stage-shuffle)
#   E: matched-random-gene-set null for COLOC + DGIdb

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
RD <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/mcp/reviewer_defense")
OUT <- FIGS_MCP_DIR

# ----- A: cophenetic across k ----------------------------------------------
stab <- fread(file.path(RD, "stability_across_k.tsv"))
p_a <- ggplot(stab, aes(k, cophenetic_corr)) +
  geom_line(linewidth = 0.5, color = masld_colors$up) +
  geom_point(size = 2, color = masld_colors$up) +
  geom_hline(yintercept = 0.95, linetype = 2, color = "gray40") +
  ylim(0.9, 1.0) +
  labs(x = "k (number of programs)", y = "Cophenetic correlation") +
  theme_minimal(base_size = 6)
ggsave(file.path(OUT, "figSmcp_a_stability_across_k.pdf"), p_a, width = 4, height = 3, device = cairo_pdf)
message("[caption] cNMF stability across k")
cat("[figS] panel A written\n")

# ----- B: entropy distribution ---------------------------------------------
ent <- fread(file.path(RD, "entropy_distribution.tsv"))
p_b <- ggplot(ent, aes(entropy_bits)) +
  geom_histogram(bins = 12, fill = masld_colors$up, alpha = 0.7) +
  geom_vline(xintercept = 1.3, linetype = 2, color = "gray40") +
  annotate("text", x = 1.35, y = 3, label = "H = 1.3\n(shared threshold)",
           size = GEOM_TEXT_6PT, hjust = 0, color = "black") +
  geom_vline(xintercept = log2(5), linetype = 3, color = "gray60") +
  annotate("text", x = log2(5) + 0.02, y = 3,
           label = "log2(5)\n= max entropy", size = GEOM_TEXT_6PT, hjust = 0, color = "black") +
  labs(x = "Cross-cell-type entropy (bits)", y = "Programs") +
  theme_minimal(base_size = 6)
ggsave(file.path(OUT, "figSmcp_b_entropy_distribution.pdf"), p_b, width = 5, height = 3, device = cairo_pdf)
message("[caption] Program entropy distribution (k=16)")
cat("[figS] panel B written\n")

# ----- C: dataset-variance vs stage-effect ---------------------------------
dc <- fread(file.path(BASE, "Analysis/SingleCell/results_gpu_v2/mcp/reviewer_defense/dataset_confound.tsv"))
dc[, program_label := sub("cnmf_global_k16_P", "P", program)]
dc[, batch_artifact := dataset_frac_var > 0.7]
dc[, stage_sig_adj := !is.na(stage_p_adj_ds) & stage_p_adj_ds < 0.05]
p_c <- ggplot(dc, aes(dataset_frac_var, -log10(pmax(stage_p_adj_ds, 1e-15)))) +
  geom_point(aes(color = stage_sig_adj, shape = batch_artifact), size = 2.5, alpha = 0.85) +
  ggrepel::geom_text_repel(aes(label = program_label), size = GEOM_TEXT_6PT, max.overlaps = 20) +
  geom_vline(xintercept = 0.7, linetype = 2, color = "gray40") +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "gray40") +
  scale_color_manual(values = c(`TRUE` = "#C2185B", `FALSE` = "gray70")) +
  scale_shape_manual(values = c(`TRUE` = 4, `FALSE` = 16)) +
  labs(x = "Fraction of variance explained by dataset",
       y = "-log10(stage p, adjusted for dataset)",
       color = "Stage sig after adj",
       shape = "Batch-dominated (>70%)") +
  theme_minimal(base_size = 6)
ggsave(file.path(OUT, "figSmcp_c_dataset_confound.pdf"), p_c, width = 6, height = 4, device = cairo_pdf)
message("[caption] Dataset confound vs stage effect per program")
cat("[figS] panel C written\n")

# ----- D: changepoint null distribution -------------------------------------
cp_null_f <- file.path(RD, "changepoint_shuffle_null_R.tsv")
if (file.exists(cp_null_f)) {
  cp <- fread(cp_null_f)
  p_d <- ggplot(cp, aes(y = program, x = null_frac_near_2)) +
    geom_col(aes(fill = obs_near_2)) +
    geom_vline(xintercept = 0.3, linetype = 2, color = "gray40") +
    scale_fill_manual(values = c(`TRUE` = masld_colors$up, `FALSE` = "gray70")) +
    labs(x = "Fraction of null fits with breakpoint near stage 2.0 (±0.1)",
         y = NULL, fill = "Observed near 2.0") +
    theme_minimal(base_size = 6) +
    theme(axis.text.y = element_text(size = 6))
  ggsave(file.path(OUT, "figSmcp_d_changepoint_null.pdf"), p_d, width = 7, height = 4, device = cairo_pdf)
  message("[caption] Changepoint null: F2 clustering is mechanical -- null fraction ~30% reveals segmented-regression bias on 4-point ordinal")
  cat("[figS] panel D written\n")
}

# ----- E: matched-null enrichment ------------------------------------------
coloc_null <- fread(file.path(RD, "coloc_enrichment_null_k16.tsv"))
drug_null <- fread(file.path(RD, "drug_enrichment_null_k16.tsv"))
coloc_null[, set := "COLOC"]; drug_null[, set := "DGIdb"]
null_combo <- rbindlist(list(coloc_null[, .(program, set, permutation_p, permutation_q_bh)],
                             drug_null[, .(program, set, permutation_p, permutation_q_bh)]))
null_combo[, program := factor(program, levels = as.character(1:20))]
p_e <- ggplot(null_combo, aes(program, -log10(pmax(permutation_p, 0.001)))) +
  geom_col(aes(fill = permutation_q_bh < 0.05)) +
  facet_wrap(~ set, ncol = 1) +
  scale_fill_manual(values = c(`TRUE` = masld_colors$up, `FALSE` = "gray70")) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "gray40") +
  labs(x = "cNMF program", y = "-log10(permutation p) [matched-random null, 1000 iter]",
       fill = "q < 0.05") +
  theme_minimal(base_size = 6)
ggsave(file.path(OUT, "figSmcp_e_matched_null_enrichment.pdf"), p_e, width = 7, height = 5, device = cairo_pdf)
message("[caption] Matched-random-gene-set null for enrichment")
cat("[figS] panel E written\n")

cat("[figS_mcp_reviewer_defense] DONE\n")
