#!/usr/bin/env Rscript
# figS04_pqtl_eqtl_crossref.R
# ─────────────────────────────────────────────────────────────────────────────
# Supplementary panel (Fig 2 / Fig 5 companion): protein-QTL vs expression-QTL
# colocalization for liver-proteogenomic candidates reported by Gobeil et al.
# 2026 (medRxiv 10.64898/2026.06.04.26354903). A TWO-EVIDENCE-AXIS scatter:
#     x = our expression-QTL colocalization (max SuSiE/ABF PP.H4, this study)
#     y = Gobeil protein-QTL colocalization (PPH4)
# Gate lines: x = 0.5 (our COLOC gate), y = 0.80 (Gobeil's pQTL prioritization).
# MTARC1 sits top-left: the source-reported protein-QTL PPH4 is 1.00 while the
# evaluated expression-QTL map has no PP.H4 > 0.5 support. GCKR is flagged as direction-
# discordant (eQTL and pQTL effect signs oppose).
#
# EXTERNAL, cited annotation only — NOT re-derived by us and NOT an input to the
# convergence/heuristic score. Reads the frozen join built by
# scripts/manuscript/build_external_pqtl_annotation.R.
#
# Conventions: PDF, 6pt, all text BLACK, gene names italic, no title/subtitle/
# in-plot annotations (caption via message()); control/neither gray #9E9E9E.
# Output: figures/main/fig2_genetics/panels/FigS2P_pqtl_eqtl_crossref.pdf (+ _source.csv)
# Env: rnaseq
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/fig2_promoted_coloc_context.R"))

ann <- fread(file.path(BASE,
  "RNA-seq/results/multi_evidence/external_pqtl/pqtl_external_annotation.tsv"))

num <- function(x) suppressWarnings(as.numeric(x))
d <- ann[!is.na(num(gobeil_pqtl_coloc_pph4)) & gene != "DUSP23"]  # 6 causal-relevant genes
d <- merge(d, load_fig2_promoted_context(BASE, d$gene), by = "gene", all.x = TRUE)
d[, x_expr := all_best]
d[, y_prot := num(gobeil_pqtl_coloc_pph4)]
d[, class := factor(all_support,
  levels = c("multi-signal COLOC", "single-signal COLOC only",
             "evaluated without PP.H4 > 0.5 support", "not evaluable"))]
d[, discordant := layer_relationship == "direction_discordant"]

class_cols <- c("multi-signal COLOC" = "#1565C0",
                "single-signal COLOC only" = "#90CAF9",
                "evaluated without PP.H4 > 0.5 support" = "#9E9E9E",
                "not evaluable" = "white")

p <- ggplot(d, aes(x_expr, y_prot)) +
  geom_hline(yintercept = 0.80, linetype = "dashed", linewidth = 0.25, colour = "grey55") +
  geom_vline(xintercept = 0.50, linetype = "dashed", linewidth = 0.25, colour = "grey55") +
  geom_point(aes(colour = class, shape = discordant), size = 1.9, stroke = 0.5) +
  ggrepel::geom_text_repel(aes(label = gene), fontface = "italic", size = GEOM_TEXT_6PT,
                           colour = "black", min.segment.length = 0,
                           segment.size = 0.2, max.overlaps = Inf, seed = 42) +
  scale_colour_manual(values = class_cols, name = "Expression-QTL support", drop = FALSE) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17),
                     labels = c(`FALSE` = "concordant sign", `TRUE` = "eQTL/pQTL sign discordant"),
                     name = NULL) +
  scale_x_continuous(limits = c(-0.02, 1.02), breaks = c(0, 0.5, 1),
                     expand = expansion(mult = 0.02)) +
  scale_y_continuous(limits = c(-0.02, 1.02), breaks = c(0, 0.5, 0.8, 1),
                     expand = expansion(mult = 0.02)) +
  labs(x = "Expression-QTL colocalization (this study, PP.H4)",
       y = "Protein-QTL colocalization (Gobeil et al. 2026, PPH4)") +
  theme_masld(base_size = 6) +
  theme(legend.key.size = unit(3, "mm"), legend.position = "right")

out_dir <- Sys.getenv("FIG2_SUPP_OUT_DIR",
  unset = file.path(BASE, "figures/main/fig2_genetics/panels"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_pdf <- file.path(out_dir, "FigS2P_pqtl_eqtl_crossref.pdf")
fwrite(d[, .(gene, expr_qtl_coloc_ours = x_expr, protein_qtl_coloc_gobeil = y_prot,
             multi_signal_pph4 = all_susie, single_signal_pph4 = all_abf,
             expression_qtl_support = as.character(class),
             direct_multi_signal_pph4 = direct_susie,
             direct_single_signal_pph4 = direct_abf,
             direct_support_state = direct_support,
             enzyme_multi_signal_pph4 = enzyme_susie,
             enzyme_single_signal_pph4 = enzyme_abf,
             enzyme_support_state = enzyme_support,
             direction_discordant = discordant)],
       file.path(out_dir, "FigS2P_pqtl_eqtl_crossref_source.csv"))
save_fig(p, out_pdf, width = 3.6, height = 3.0)

message("Wrote ", out_pdf)
message("CAPTION: Source-reported protein-QTL colocalization (Gobeil et al. 2026) ",
        "versus the promoted expression-QTL COLOC portfolio. Color separates multi-signal, ",
        "single-signal-only, evaluated unsupported, and not-evaluable expression-QTL states; ",
        "shape retains the source-reported eQTL/pQTL direction comparison. External annotation; ",
        "not a Resource evidence-class input.")
