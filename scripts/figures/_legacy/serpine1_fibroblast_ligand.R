#!/usr/bin/env Rscript
# ==============================================================================
# POSITIVE CONTROL: SERPINE1 (PAI-1), an activated-HSC fibroblast-derived ligand.
#
# KEY MESSAGE
#   "SERPINE1, a top fibroblast-derived ligand by LIANA (0.916), is upregulated
#    in disease and gains spatial organization — recovering the activated-HSC
#    PAI-1 axis as a method positive control."
#
# HONESTY GUARD-RAILS (baked into the legend via message() to stdout):
#   * SERPINE1 is "a top" fibroblast ligand, NOT the maximum: by LIANA score it
#     ranks 3rd of 10 fibroblast-sourced ligands (LAMA4 0.972, NUCB2 0.968 rank
#     above). What distinguishes it is the COMBINATION of high ligand score AND
#     high bulk-disease upregulation (log2FC +0.54) — LAMA4/NUCB2 are near-flat.
#     The landscape shows that combination honestly rather than asserting #1.
#   * Spatial Moran's I 0.059 -> 0.541 is INCREASED ORGANIZATION only, NOT a
#     disease-direction effect (GSE192741, n=5 donors, batch-flagged).
#   * Endothelium is a biological co-source of PAI-1 (Kim 2024 J Hepatol).
#   * COLOC PP.H4 = 0.14 (abf, UKBB_GGT; no SuSiE) — NO genetic pillar.
#   * Causality for HSC-derived PAI-1 was established by perturbation in
#     Kim 2024 J Hepatol (PMID 39522884); WE provide ASSOCIATION not
#     perturbation, hence "positive control" not novel causal claim.
#
# Panels:
#   i.  Fibroblast-derived ligands: bulk log2FC (x) x LIANA score (y); SERPINE1 hero
#   ii. SERPINE1 spatial autocorrelation (Moran's I), healthy -> steatotic
#
# Bulk numbers are the canonical limma-voom quality-weighted C2 channel
# (secretome_chain producer file `bulk_lfc`), NOT the retired dream method.
#
# Output: figures/main/fig4_validation/fig4e_serpine1_fibroblast_ligand.pdf
# Env:    rnaseq
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(ggrepel)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# ── Data: fibroblast-ligand landscape (LIANA producer file) ──────────────────
# Producer file (NOT an inverted consumer): max_score_diff = LIANA ligand score,
# bulk_lfc = canonical C2 disease-vs-control log2FC, n_masld/n_ctrl_enriched =
# MASLD-vs-control pair counts. Ledger rows 39 (LIANA 0.916) + 20 (bulk +0.522;
# secretome file uses a separately-computed +0.544 — both VERIFIED).
secretome <- read.csv(file.path(BASE,
  "RNA-seq/results/secretome_chain/secretome_chain_triple_concordant_MASLD_up.csv"),
  stringsAsFactors = FALSE) %>%
  mutate(gene = ligand_complex)

# One row per gene (best LIANA score), restricted to fibroblast-sourced ligands.
fib <- secretome %>%
  mutate(src = ifelse(!is.na(primary_source_ct) & primary_source_ct != "",
                      primary_source_ct, a1_primary)) %>%
  filter(grepl("fibro", src, ignore.case = TRUE),
         !is.na(max_score_diff), !is.na(bulk_lfc)) %>%
  group_by(gene) %>% slice_max(max_score_diff, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(is_focal = gene == "SERPINE1")

serp <- fib %>% filter(is_focal)
stopifnot(nrow(serp) == 1)
serp_rank <- which(fib$gene[order(-fib$max_score_diff)] == "SERPINE1")

# ── Data: per-condition spatial autocorrelation (single source of truth) ─────
# Ledger row 32 / 49: GSE192741 per-condition Moran's I. The atlas single-value
# spatial col (0.242) is method-inconsistent — do NOT cite it.
autocorr_h <- read.csv(file.path(BASE,
  "Analysis/Spatial/results/coexpression/spatial_autocorr_Healthy.csv"),
  stringsAsFactors = FALSE)
autocorr_s <- read.csv(file.path(BASE,
  "Analysis/Spatial/results/coexpression/spatial_autocorr_Steatotic.csv"),
  stringsAsFactors = FALSE)
serp_moran_h <- autocorr_h$C[autocorr_h$Gene == "SERPINE1"][1]   # 0.059 healthy
serp_moran_s <- autocorr_s$C[autocorr_s$Gene == "SERPINE1"][1]   # 0.541 steatotic

# ── Data: COLOC (honesty — no genetic pillar) ────────────────────────────────
# Ledger row 11: SERPINE1 abf PP.H4 = 0.138 (UKBB_GGT), no SuSiE convergence.
coloc_all <- read.csv(file.path(BASE,
  "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv"),
  stringsAsFactors = FALSE)
serp_coloc <- max(coloc_all$PP.H4.abf[coloc_all$gene == "SERPINE1"], na.rm = TRUE)  # 0.138

# ══════════════════════════════════════════════════════════════════════════════
# Panel i — fibroblast-derived ligand landscape
# ══════════════════════════════════════════════════════════════════════════════
# SERPINE1 occupies the upper-RIGHT (high ligand score AND high disease LFC).
# LAMA4 / NUCB2 sit higher on y but near-flat on x. Direct label, no legend.
p_land <- ggplot(fib, aes(x = bulk_lfc, y = max_score_diff)) +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "grey75") +
  geom_point(aes(color = is_focal, size = is_focal)) +
  geom_text_repel(aes(label = gene,
                      fontface = ifelse(is_focal, "bold.italic", "italic")),
                  color = "black",
                  size = PUB_GEOM_TEXT, box.padding = 0.30, segment.size = 0.2,
                  segment.color = "grey70", max.overlaps = Inf, seed = 1,
                  show.legend = FALSE) +
  scale_color_manual(values = c(`TRUE` = masld_colors[["up"]], `FALSE` = "#C9A0B4"),
                     guide = "none") +
  scale_size_manual(values = c(`TRUE` = 3.0, `FALSE` = 1.4), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.10, 0.18))) +
  scale_y_continuous(limits = c(0.70, 1.0)) +
  labs(x = "Bulk mRNA log2FC (disease vs control)",
       y = "LIANA ligand score") +
  theme_masld() + theme_pub()

# ══════════════════════════════════════════════════════════════════════════════
# Panel ii — SERPINE1 spatial autocorrelation (organization), healthy -> steatotic
# ══════════════════════════════════════════════════════════════════════════════
spat_df <- data.frame(
  cond  = factor(c("Healthy", "Steatotic"), levels = c("Healthy", "Steatotic")),
  moran = c(serp_moran_h, serp_moran_s),
  fill  = c(masld_colors[["control"]], masld_colors[["up"]])
)
p_spat <- ggplot(spat_df, aes(x = cond, y = moran, fill = fill)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = sprintf("%.2f", moran)),
            vjust = -0.4, size = PUB_GEOM_TEXT, color = "grey20") +
  scale_fill_identity() +
  scale_y_continuous(limits = c(0, max(spat_df$moran) * 1.14),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Spatial autocorrelation\n(Moran's I)") +
  theme_masld() + theme_pub() +
  theme(axis.text.x = element_text(size = PUB_AXIS_TEXT + 0.5))

# ══════════════════════════════════════════════════════════════════════════════
# Legend / provenance — emitted to stdout (stats live in the legend, NOT panel)
# ══════════════════════════════════════════════════════════════════════════════
message(strrep("=", 78))
message("FIGURE LEGEND (paste into manuscript; stats live here, not on the panel)")
message(strrep("=", 78))
message(sprintf(
"SERPINE1 (PAI-1) as a method positive control. (i) Among fibroblast-derived
secreted ligands, SERPINE1 is a TOP ligand by LIANA score (%.3f; rank %d of %d
fibroblast ligands — LAMA4 0.972 and NUCB2 0.968 rank above it) and is the
fibroblast ligand that pairs a high ligand score with strong bulk-disease
upregulation (log2FC %+.2f, padj %.1e; canonical limma-voom-qw C2; MASLD-enriched
in %d vs control %d of %d LIANA pairs). (ii) Its spatial structure strengthens
from Healthy to Steatotic liver (Moran's I %.3f -> %.3f; GSE192741).",
  serp$max_score_diff[1], serp_rank, nrow(fib),
  serp$bulk_lfc[1], serp$bulk_padj[1],
  serp$n_masld_enriched[1], serp$n_ctrl_enriched[1], serp$n_pairs[1],
  serp_moran_h, serp_moran_s))
message("")
message("CAVEATS (mandatory):")
message(sprintf(
" - Moran's I %.3f -> %.3f reflects INCREASED SPATIAL ORGANIZATION only, NOT a
   disease-direction effect; GSE192741 has n=5 donors (2 Healthy / 3 Steatotic
   JBO) and spot-level autocorrelation makes per-spot permutation anticonservative
   (batch-flagged).", serp_moran_h, serp_moran_s))
message(" - SERPINE1 is 'a top' fibroblast ligand, not the maximum (LAMA4 0.972,
   NUCB2 0.968 score higher but are near-flat in bulk).")
message(" - Endothelium is a biological CO-SOURCE of PAI-1; this is not a
   fibroblast-exclusive ligand.")
message(sprintf(
" - COLOC PP.H4 = %.2f (abf, UKBB_GGT; no SuSiE convergence) — NO genetic pillar
   for SERPINE1. The SERPINE1 story is cell-cell-communication + spatial, not
   colocalization.", serp_coloc))
message(" - Causality for activated-HSC-derived PAI-1 was established by
   perturbation in Kim et al. 2024 J Hepatol (PMID 39522884). WE provide
   association (LIANA + bulk + spatial), not perturbation — hence this panel is
   a METHOD POSITIVE CONTROL recovering a known axis, not a novel causal claim.")
message(strrep("=", 78))

# ══════════════════════════════════════════════════════════════════════════════
# Assemble + save
# ══════════════════════════════════════════════════════════════════════════════
p_out <- (p_land | p_spat) + plot_layout(widths = c(2.2, 1))

out <- file.path(FIG4_DIR, "fig4e_serpine1_fibroblast_ligand.pdf")
cairo_pdf(out, width = fig_full_width * 0.92, height = 2.5, onefile = FALSE)
print(p_out)
invisible(dev.off())
message("Saved: ", out)
