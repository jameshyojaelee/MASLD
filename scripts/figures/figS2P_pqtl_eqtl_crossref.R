#!/usr/bin/env Rscript
# figS04_pqtl_eqtl_crossref.R
# ─────────────────────────────────────────────────────────────────────────────
# Supplementary panel (Fig 2 / Fig 5 companion): protein-QTL vs expression-QTL
# colocalization for the liver-proteogenomic causal effectors of Gobeil et al.
# 2026 (medRxiv 10.64898/2026.06.04.26354903). A TWO-EVIDENCE-AXIS scatter:
#     x = our expression-QTL colocalization (max SuSiE/ABF PP.H4, this study)
#     y = Gobeil protein-QTL colocalization (PPH4)
# Gate lines: x = 0.5 (our COLOC gate), y = 0.80 (Gobeil's pQTL prioritization).
# The point of the panel: MTARC1 sits top-left — protein-QTL colocalizes (PPH4=1.00)
# where expression-QTL is null (0.045) — the complementary causal layer this
# expression-based resource does not generate. GCKR is flagged as direction-
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

ann <- fread(file.path(BASE,
  "RNA-seq/results/multi_evidence/external_pqtl/pqtl_external_annotation.tsv"))

num <- function(x) suppressWarnings(as.numeric(x))
d <- ann[!is.na(num(gobeil_pqtl_coloc_pph4)) & gene != "DUSP23"]  # 6 causal-relevant genes
d[, x_expr := pmax(fcoalesce(num(our_max_susie_pp4), 0),
                   fcoalesce(num(our_max_abf_pp4), 0))]
d[, y_prot := num(gobeil_pqtl_coloc_pph4)]
d[, class := factor(fifelse(our_primary_class == "genetic_only", "Genetic-only",
                    fifelse(our_primary_class == "disease_state_only", "Disease-state-only",
                            "Neither")),
                    levels = c("Genetic-only", "Disease-state-only", "Neither"))]
d[, discordant := layer_relationship == "direction_discordant"]

class_cols <- c("Genetic-only" = "#c0508d",        # Liang magenta
                "Disease-state-only" = "#3a7ec2",   # Liang blue
                "Neither" = "#9E9E9E")              # locked control/neutral gray

p <- ggplot(d, aes(x_expr, y_prot)) +
  geom_hline(yintercept = 0.80, linetype = "dashed", linewidth = 0.25, colour = "grey55") +
  geom_vline(xintercept = 0.50, linetype = "dashed", linewidth = 0.25, colour = "grey55") +
  geom_point(aes(colour = class, shape = discordant), size = 1.9, stroke = 0.5) +
  ggrepel::geom_text_repel(aes(label = gene), fontface = "italic", size = 2.1,
                           colour = "black", min.segment.length = 0,
                           segment.size = 0.2, max.overlaps = Inf, seed = 42) +
  scale_colour_manual(values = class_cols, name = "Our class", drop = FALSE) +
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

out_dir <- file.path(BASE, "figures/main/fig2_genetics/panels")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_pdf <- file.path(out_dir, "FigS2P_pqtl_eqtl_crossref.pdf")
fwrite(d[, .(gene, expr_qtl_coloc_ours = x_expr, protein_qtl_coloc_gobeil = y_prot,
             our_primary_class = as.character(class), direction_discordant = discordant)],
       file.path(out_dir, "FigS2P_pqtl_eqtl_crossref_source.csv"))
save_fig(p, out_pdf, width = 3.6, height = 3.0)

message("Wrote ", out_pdf)
message("CAPTION: Protein-QTL (Gobeil et al. 2026, liver mass-spec pQTL) vs our ",
        "expression-QTL colocalization for their liver-proteogenomic causal effectors. ",
        "MTARC1 colocalizes at the protein level (PPH4=1.00) where our expression map is null ",
        "(PP.H4=0.045, Neither class); GCKR colocalizes on expression but with a protein-level ",
        "effect direction opposite the transcript-level one. External cited annotation; not a ",
        "convergence-score input.")