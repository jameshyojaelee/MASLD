#!/usr/bin/env Rscript
# ==============================================================================
# FIBROBLAST-DERIVED LIGAND validation example: THBS2 (thrombospondin-2).
#
# Replaces the retired SERPINE1 panel: under the CORRECTED CCC sign convention
# (2026-06-20) SERPINE1's top MASLD-enriched sender reassigned away from
# fibroblasts (an unreliable pooled-LIANA within-condition-normalised source
# call). THBS2 is robustly fibroblast-sourced in the corrected secretome and a
# stronger, cleaner cross-modal example.
#
# KEY MESSAGE
#   "THBS2, a fibroblast-derived secreted ligand, pairs the STRONGEST bulk
#    disease upregulation among top fibroblast ligands (log2FC +1.14, padj
#    3.6e-35) with a high LIANA cell-cell-communication score and plasma (Olink)
#    detectability — recovering the known THBS2 NASH-fibrosis serum-biomarker
#    axis as a method positive control."
#
# HONESTY GUARD-RAILS (in the stdout legend, not on the panel):
#   * THBS2 is the fibroblast ligand with the strongest bulk-disease upregulation
#     among the high-LIANA set; by LIANA score alone HGF/VEGFC rank above it (but
#     with near-flat bulk). The landscape shows the combination honestly.
#   * Plasma evidence = Olink-panel presence + established clinical NASH-fibrosis
#     biomarker (literature); our DIA-MS proteomics did not separately quantify it.
#   * Fibroblasts / activated HSC are a major but not exclusive source.
#   * COLOC reported honestly — the THBS2 story is CCC + bulk + plasma, not
#     colocalization.
#
# Panel: fibroblast-derived ligands, bulk log2FC (x) x LIANA ligand score (y);
#        THBS2 the hero. Bulk = canonical limma-voom-qw C2 (`bulk_lfc`).
# Output: figures/main/fig4_validation/fig4e_thbs2_fibroblast_ligand.pdf
# Env:    rnaseq
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

FOCAL <- "THBS2"

# ── Data: fibroblast-ligand landscape (LIANA secretome producer file) ─────────
# Producer file (NOT an inverted consumer): max_score_diff = LIANA ligand score
# (corrected sign, score_diff>0 == MASLD-enriched), bulk_lfc = canonical C2
# disease-vs-control log2FC.
secretome <- read.csv(file.path(BASE,
  "RNA-seq/results/secretome_chain/secretome_chain_triple_concordant_MASLD_up.csv"),
  stringsAsFactors = FALSE) %>%
  mutate(gene = ligand_complex)

# One row per gene (best LIANA score), restricted to fibroblast-sourced ligands.
fib <- secretome %>%
  mutate(src = ifelse(!is.na(primary_source_ct) & primary_source_ct != "",
                      as.character(primary_source_ct), as.character(a1_primary))) %>%
  filter(grepl("fibro", src, ignore.case = TRUE),
         !is.na(max_score_diff), !is.na(bulk_lfc)) %>%
  group_by(gene) %>% slice_max(max_score_diff, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(is_focal = gene == FOCAL)

foc <- fib %>% filter(is_focal)
stopifnot(nrow(foc) == 1)
foc_rank <- which(fib$gene[order(-fib$max_score_diff)] == FOCAL)

# ── COLOC (honesty — report whatever it is) ──────────────────────────────────
coloc_all <- read.csv(file.path(BASE,
  "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv"),
  stringsAsFactors = FALSE)
foc_coloc <- suppressWarnings(max(coloc_all$PP.H4.abf[coloc_all$gene == FOCAL], na.rm = TRUE))
if (!is.finite(foc_coloc)) foc_coloc <- NA_real_

# ══════════════════════════════════════════════════════════════════════════════
# Panel — fibroblast-derived ligand landscape (THBS2 the upper-right hero)
# ══════════════════════════════════════════════════════════════════════════════
y_lo <- floor(min(fib$max_score_diff, na.rm = TRUE) * 20) / 20
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
  scale_y_continuous(limits = c(y_lo, 1.0)) +
  labs(x = "Bulk mRNA log2FC (disease vs control)",
       y = "LIANA ligand score") +
  theme_masld() + theme_pub()

# ══════════════════════════════════════════════════════════════════════════════
# Legend / provenance — emitted to stdout (stats live here, NOT on the panel)
# ══════════════════════════════════════════════════════════════════════════════
message(strrep("=", 78))
message("FIGURE LEGEND (paste into manuscript; stats live here, not on the panel)")
message(strrep("=", 78))
message(sprintf(
"THBS2 (thrombospondin-2) as a fibroblast-derived-ligand positive control. Among
fibroblast-sourced secreted ligands, THBS2 shows the STRONGEST bulk-disease
upregulation (log2FC %+.2f, padj %.1e; canonical limma-voom-qw C2; MASLD-enriched
in %d vs control %d of %d LIANA pairs) paired with a high LIANA cell-cell-
communication ligand score (%.3f; rank %d of %d fibroblast ligands by LIANA score).
THBS2 is a recognised circulating NASH-fibrosis biomarker and is present in the
Olink plasma panel — recovering a known fibrosis axis as a method positive control.",
  foc$bulk_lfc[1], foc$bulk_padj[1],
  foc$n_masld_enriched[1], foc$n_ctrl_enriched[1], foc$n_pairs[1],
  foc$max_score_diff[1], foc_rank, nrow(fib)))
message("")
message("CAVEATS (mandatory):")
message(" - THBS2 is the fibroblast ligand with the strongest bulk upregulation, NOT the\n   maximum LIANA score (HGF 0.974 / VEGFC 0.943 score higher on LIANA but are\n   near-flat in bulk). The landscape shows the combination honestly.")
message(" - Plasma evidence = Olink-panel presence + established clinical NASH-fibrosis\n   biomarker (literature); our DIA-MS proteomics did not separately quantify THBS2.")
message(" - Fibroblasts / activated HSC are a major but not exclusive source of THBS2.")
message(sprintf(
" - COLOC PP.H4 (abf) = %s — the THBS2 story is cell-cell-communication + bulk +\n   plasma, NOT colocalization (no genetic pillar claimed).",
  ifelse(is.na(foc_coloc), "NA (no colocalization signal)", sprintf("%.2f", foc_coloc))))
message(strrep("=", 78))

# ══════════════════════════════════════════════════════════════════════════════
# Save
# ══════════════════════════════════════════════════════════════════════════════
out <- file.path(FIG4_DIR, "fig4e_thbs2_fibroblast_ligand.pdf")
ggsave(out, p_land, width = 4.8, height = 3.5, useDingbats = FALSE)
message("Wrote: ", out)
