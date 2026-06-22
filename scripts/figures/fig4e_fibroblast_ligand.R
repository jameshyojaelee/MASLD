#!/usr/bin/env Rscript
# ==============================================================================
# FIBROBLAST-DERIVED LIGAND validation example for Fig 4e (parameterised hero).
#
# Replaces the retired SERPINE1 panel: under the CORRECTED CCC sign convention
# (2026-06-20) SERPINE1's top MASLD-enriched sender reassigned away from
# fibroblasts (an unreliable pooled-LIANA within-condition-normalised source
# call). This script highlights any robustly-fibroblast-sourced ligand as the
# hero so a candidate can be picked for the final figure.
#
# Choose the hero with the FIG4E_FOCAL env var (default THBS2). Candidates worth
# comparing (all in the triple-concordant fibroblast set, all Olink-present):
#   THBS2  bulk +1.14 / LIANA 0.848  — strongest bulk; NASH-fibrosis biomarker
#   EFEMP1 bulk +1.08 / LIANA 0.836  — bulk+CCC+measured DIA-MS plasma protein
#   COL4A1 bulk +0.56 / LIANA 0.748  — basement-membrane collagen IV
#   HGF    bulk +0.25 / LIANA 0.974  — highest CCC, weak bulk (alt profile)
#
# Panel: fibroblast-derived ligands, bulk log2FC (x) x LIANA ligand score (y);
#        FOCAL the hero. Bulk = canonical limma-voom-qw C2 (`bulk_lfc`).
# Output: figures/main/fig4_validation/fig4e_fibroblast_ligand_<FOCAL>.pdf
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

FOCAL <- Sys.getenv("FIG4E_FOCAL", "THBS2")

# Per-gene biomarker context for the legend (generic fallback otherwise).
BIO_CONTEXT <- list(
  THBS2  = "a recognised circulating NASH-fibrosis biomarker",
  EFEMP1 = "fibulin-3, an established liver-fibrosis matricellular marker (and the one fibroblast ligand here with a measured DIA-MS plasma elevation)",
  COL4A1 = "basement-membrane collagen IV, a core fibrosis ECM component",
  HGF    = "hepatocyte growth factor, a canonical mesenchymal->hepatocyte signal",
  VEGFC  = "a lymph/endothelial growth factor",
  COL4A2 = "basement-membrane collagen IV",
  SLIT2  = "an axon-guidance ligand re-deployed in fibrosis"
)
bio <- if (!is.null(BIO_CONTEXT[[FOCAL]])) BIO_CONTEXT[[FOCAL]] else "a fibroblast-derived secreted ligand"

# ── Data: fibroblast-ligand landscape (LIANA secretome producer file) ─────────
secretome <- read.csv(file.path(BASE,
  "RNA-seq/results/secretome_chain/secretome_chain_triple_concordant_MASLD_up.csv"),
  stringsAsFactors = FALSE) %>%
  mutate(gene = ligand_complex)

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
foc_rank_liana <- which(fib$gene[order(-fib$max_score_diff)] == FOCAL)
foc_rank_bulk  <- which(fib$gene[order(-fib$bulk_lfc)] == FOCAL)
is_top_bulk    <- foc_rank_bulk == 1L

coloc_all <- read.csv(file.path(BASE,
  "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv"),
  stringsAsFactors = FALSE)
foc_coloc <- suppressWarnings(max(coloc_all$PP.H4.abf[coloc_all$gene == FOCAL], na.rm = TRUE))
if (!is.finite(foc_coloc)) foc_coloc <- NA_real_

# ══════════════════════════════════════════════════════════════════════════════
# Panel — fibroblast-derived ligand landscape (FOCAL highlighted)
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

# ── Legend / provenance to stdout (stats live here, NOT on the panel) ─────────
bulk_phrase <- if (is_top_bulk) {
  sprintf("the STRONGEST bulk-disease upregulation among fibroblast ligands (log2FC %+.2f, padj %.1e)", foc$bulk_lfc[1], foc$bulk_padj[1])
} else {
  sprintf("strong bulk-disease upregulation (log2FC %+.2f, padj %.1e; rank %d of %d by bulk)", foc$bulk_lfc[1], foc$bulk_padj[1], foc_rank_bulk, nrow(fib))
}
message(strrep("=", 78))
message(sprintf("FIGURE LEGEND — hero = %s", FOCAL))
message(strrep("=", 78))
message(sprintf(
"%s (%s) as a fibroblast-derived-ligand positive control. Among fibroblast-sourced
secreted ligands, %s shows %s, paired with a LIANA cell-cell-communication ligand
score of %.3f (rank %d of %d fibroblast ligands by LIANA; canonical limma-voom-qw C2;
MASLD-enriched in %d vs control %d of %d LIANA pairs). %s is present in the Olink
plasma panel.",
  FOCAL, bio, FOCAL, bulk_phrase, foc$max_score_diff[1], foc_rank_liana, nrow(fib),
  foc$n_masld_enriched[1], foc$n_ctrl_enriched[1], foc$n_pairs[1], FOCAL))
message("")
message("CAVEATS (mandatory):")
if (!is_top_bulk)
  message(sprintf(" - %s scores high on LIANA but is NOT the strongest in bulk; the powered\n   disease signal (bulk/plasma) is the load-bearing axis.", FOCAL))
message(" - LIANA source assignment is from the pooled (within-condition-normalised)\n   CCC layer; treat the cell-type SOURCE as supporting, not definitive.")
message(" - Fibroblasts / activated HSC are a major but not exclusive source.")
message(sprintf(
" - COLOC PP.H4 (abf) = %s — the story is cell-cell-communication + bulk + plasma,\n   NOT colocalization (no genetic pillar claimed).",
  ifelse(is.na(foc_coloc), "NA (no colocalization signal)", sprintf("%.2f", foc_coloc))))
message(strrep("=", 78))

out <- file.path(FIG4_DIR, sprintf("fig4e_fibroblast_ligand_%s.pdf", FOCAL))
ggsave(out, p_land, width = 4.8, height = 3.5, useDingbats = FALSE)
message("Wrote: ", out)
