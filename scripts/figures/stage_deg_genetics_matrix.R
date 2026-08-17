#!/usr/bin/env Rscript
# KEY MESSAGE: Stage-associated expression and exact shared-signal genetics provide distinct, traceable evidence for the same genes.

# Figure 3F. Adjacent-stage effects joined to an exact promoted COLOC signal
# pair. PP.H4 and PIP remain separate lanes; unresolved PIP is never zero-filled.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

fail <- function(...) stop(..., call. = FALSE)
sha256_file <- function(path) {
  x <- system2("sha256sum", normalizePath(path, mustWork = TRUE), stdout = TRUE)
  strsplit(x[[1L]], "[[:space:]]+")[[1L]][[1L]]
}
base <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
candidate_root <- Sys.getenv("FIGURE_CANDIDATE_ROOT", "")
stage_root <- normalizePath(Sys.getenv("STAGE_RELEASE_ROOT", ""), mustWork = TRUE)
genetics_root <- normalizePath(Sys.getenv("GENETICS_RELEASE_ROOT", ""), mustWork = TRUE)
if (!nzchar(candidate_root) || !dir.exists(candidate_root)) fail("Active FIGURE_CANDIDATE_ROOT required")
source(file.path(base, "scripts/figures/publication_theme.R"))

promotion_json <- file.path(genetics_root, "PROMOTED.json")
manifest_path <- file.path(genetics_root, "output_manifest.tsv")
pair_path <- file.path(genetics_root, "figure_signal_pairs.tsv")
if (!all(file.exists(c(promotion_json, manifest_path, pair_path)))) {
  fail("Genetics gate closed: promoted terminal manifest and figure_signal_pairs.tsv are required")
}
promotion <- paste(readLines(promotion_json, warn = FALSE), collapse = "")
if (!grepl('"status"[[:space:]]*:[[:space:]]*"promoted"', promotion) ||
    !grepl('"terminal"[[:space:]]*:[[:space:]]*true', promotion) ||
    !grepl('"canonical_promotion_authorized"[[:space:]]*:[[:space:]]*true', promotion)) {
  fail("Genetics gate closed: corrected COLOC release is not terminal and promoted")
}
manifest <- fread(manifest_path, colClasses = "character")
required_manifest <- c("relative_path", "size_bytes", "sha256")
if (!all(required_manifest %in% names(manifest))) fail("Genetics output manifest schema drift")
pair_row <- manifest[relative_path == "figure_signal_pairs.tsv"]
if (nrow(pair_row) != 1L || pair_row$sha256 != sha256_file(pair_path) ||
    as.numeric(pair_row$size_bytes) != file.info(pair_path)$size) {
  fail("Genetics signal-pair checksum mismatch")
}

required_pairs <- c(
  "gene_id_versioned", "gene_symbol", "trait_scope", "trait", "gwas_name",
  "coloc_pair_id", "PP.H4.susie", "primary_colocalized", "credible_set_id",
  "lead_variant", "pip", "linkage_status", "pooled_logFC", "pooled_FDR",
  "named_in_final_manuscript", "finalized_figure2_gene"
)
pairs <- fread(pair_path)
if (!all(required_pairs %in% names(pairs))) fail("Promoted signal-pair export schema is incomplete")
if (anyDuplicated(pairs$coloc_pair_id) || any(pairs$PP.H4.susie < 0 | pairs$PP.H4.susie > 1, na.rm = TRUE) ||
    any(pairs$pip < 0 | pairs$pip > 1, na.rm = TRUE)) fail("Invalid signal-pair values")

stage_path <- file.path(candidate_root, "analysis", "stage_extensions", "stage_all_gene_results.tsv")
stage <- fread(stage_path)
transitions <- c("F0_to_F1", "F1_to_F2", "F2_to_F3", "F3_to_F4")
if (any(stage[, uniqueN(gene_id_versioned), by = contrast]$V1 != unique(stage$bh_family_size))) fail("Incomplete stage BH family")
eligible_ids <- unique(stage[contrast %in% transitions & FDR < 0.05, gene_id_versioned])
primary <- pairs[primary_colocalized == TRUE & gene_id_versioned %in% eligible_ids]
if (!nrow(primary)) fail("No promoted primary colocalized stage-associated genes")
setorder(primary, gene_id_versioned, -PP.H4.susie, trait_scope, trait, coloc_pair_id)
primary <- primary[, .SD[1L], by = gene_id_versioned]

matrix <- merge(
  stage[gene_id_versioned %in% primary$gene_id_versioned & contrast %in% transitions,
    .(gene_id_versioned, gene_name, transition = contrast, logFC, FDR, n_samples, n_cohorts)],
  primary, by = "gene_id_versioned", all.x = TRUE, sort = FALSE
)
if (matrix[, uniqueN(transition), by = gene_id_versioned][, any(V1 != 4L)]) fail("Incomplete four-contrast matrix")
rank_table <- matrix[, .(
  max_abs_effect = max(abs(logFC)),
  max_effect_direction = sign(logFC[which.max(abs(logFC))]),
  min_stage_fdr = min(FDR),
  max_transition = transition[which.max(abs(logFC))],
  PP.H4.susie = PP.H4.susie[[1L]],
  named_in_final_manuscript = named_in_final_manuscript[[1L]]
), by = .(gene_id_versioned, gene_symbol)]
rank_table[, transition_rank := match(max_transition, transitions)]
setorder(rank_table, transition_rank, -max_effect_direction, min_stage_fdr, gene_symbol)
rank_table[, row_order := .I]
label_ids <- rank_table[named_in_final_manuscript == TRUE][order(-max_abs_effect, -PP.H4.susie, gene_symbol)][1:min(8L, .N), gene_id_versioned]
rank_table[, row_label := fifelse(gene_id_versioned %in% label_ids, gene_symbol, "")]
matrix <- merge(matrix, rank_table[, .(gene_id_versioned, row_order, row_label)], by = "gene_id_versioned")
matrix[, transition := factor(transition, levels = transitions,
  labels = c("F0→F1", "F1→F2", "F2→F3", "F3→F4"))]
matrix[, gene_row := factor(row_order, levels = rev(rank_table$row_order))]
matrix[, display_effect := fifelse(FDR < 0.05, logFC, 0.35 * logFC)]
matrix[, significant := FDR < 0.05]

lim <- quantile(abs(matrix$logFC), 0.98, na.rm = TRUE)
heat <- ggplot(matrix, aes(transition, gene_row, fill = display_effect, alpha = significant)) +
  geom_tile(color = "white", linewidth = 0.12) +
  scale_fill_gradient2(low = masld_colors$down, mid = "white", high = masld_colors$mash,
                       midpoint = 0, limits = c(-lim, lim), oob = scales::squish,
                       name = expression(log[2]*"FC")) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.42), guide = "none") +
  scale_y_discrete(labels = setNames(rank_table$row_label, rank_table$row_order)) +
  labs(x = NULL, y = NULL) + theme_masld_compact() +
  theme(axis.text.y = element_text(size = 6, face = "italic"),
        axis.text.x = element_text(size = 6, face = "plain"), legend.position = "bottom")

lane <- unique(matrix[, .(gene_id_versioned, row_order, PP.H4.susie, pip, linkage_status, credible_set_id, lead_variant)])
lane[, gene_row := factor(row_order, levels = rev(rank_table$row_order))]
lane[, pip_display := fifelse(linkage_status == "exact_gene_trait_credible_set", pip, NA_real_)]
pp <- ggplot(lane, aes("PP.H4", gene_row, fill = PP.H4.susie)) +
  geom_tile(color = "white", linewidth = 0.12) +
  scale_fill_gradient(low = "white", high = cat_palette[[3]], limits = c(0, 1), name = "PP.H4") +
  labs(x = NULL, y = NULL) + theme_masld_compact() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "bottom")
pip_plot <- ggplot(lane, aes("PIP", gene_row, fill = pip_display)) +
  geom_tile(color = "white", linewidth = 0.12) +
  geom_point(data = lane[linkage_status != "exact_gene_trait_credible_set" | is.na(pip)],
             shape = 4, size = 0.7, stroke = 0.3, color = "grey55") +
  scale_fill_gradient(low = "white", high = cat_palette[[5]], limits = c(0, 1), na.value = "grey92", name = "PIP") +
  labs(x = NULL, y = NULL) + theme_masld_compact() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "bottom")

panel <- (heat | pp | pip_plot) + plot_layout(widths = c(4.3, 0.7, 0.7), guides = "keep")
out_dir <- file.path(candidate_root, "figure3", "panels")
supp_dir <- file.path(candidate_root, "supplementary", "figureS3", "panels")
src_dir <- file.path(candidate_root, "source_tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(src_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(out_dir, "fig3f_stage_deg_genetics_matrix.pdf"), panel,
       width = 4.45, height = 5.2, device = cairo_pdf)

scatter <- pairs[primary_colocalized == TRUE & is.finite(pooled_logFC) & is.finite(`PP.H4.susie`)]
setorder(scatter, gene_id_versioned, -`PP.H4.susie`, trait_scope, trait, coloc_pair_id)
scatter <- scatter[, .SD[1L], by = gene_id_versioned]
p_scatter <- ggplot(scatter, aes(pooled_logFC, PP.H4.susie)) +
  geom_hline(yintercept = 0.5, color = "grey70", linewidth = 0.25, linetype = 2) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.25) +
  geom_point(aes(color = pooled_FDR < 0.05), size = 0.7, alpha = 0.7) +
  scale_color_manual(values = c(`FALSE` = "#9E9E9E", `TRUE` = cat_palette[[3]]), guide = "none") +
  labs(x = expression("Pooled RNA-seq "*log[2]*"FC"), y = "COLOC PP.H4") + theme_masld_compact()
ggsave(file.path(supp_dir, "figs3_rnaseq_logfc_vs_pph4.pdf"), p_scatter,
       width = 2.5, height = 2.25, device = cairo_pdf)

pip_focus <- pairs[primary_colocalized == TRUE & linkage_status == "exact_gene_trait_credible_set" &
                     is.finite(pip) & is.finite(`PP.H4.susie`)]
p_pip <- ggplot(pip_focus,
                aes(PP.H4.susie, pip)) +
  geom_point(size = 0.8, color = cat_palette[[5]]) +
  labs(x = "COLOC PP.H4", y = "Lead GWAS SuSiE PIP") + theme_masld_compact()
ggsave(file.path(supp_dir, "figs3_exact_linked_pip.pdf"), p_pip,
       width = 2.5, height = 2.25, device = cairo_pdf)

fwrite(matrix, file.path(src_dir, "fig3f_stage_genetics_matrix.tsv"), sep = "\t")
fwrite(lane, file.path(src_dir, "fig3f_exact_pair_lanes.tsv"), sep = "\t")
fwrite(scatter, file.path(src_dir, "figs3_rnaseq_logfc_vs_pph4.tsv"), sep = "\t")
fwrite(pip_focus, file.path(src_dir, "figs3_exact_linked_pip.tsv"), sep = "\t")
cat("[saved] Figure 3F and supplementary alternatives\n")
