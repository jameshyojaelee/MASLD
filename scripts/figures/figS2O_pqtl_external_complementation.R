#!/usr/bin/env Rscript
# figS2O_pqtl_external_complementation.R
# ─────────────────────────────────────────────────────────────────────────────
# Fig 2 supplement (FigS2O): external liver protein-QTL complementation of our
# expression-QTL genetic map. Compact heatmap-table for the liver-proteogenomic
# prioritized effectors reported by Gobeil et al. 2026
# (medRxiv 10.64898/2026.06.04.26354903):
#   col 1 = our expression-QTL colocalization (max SuSiE/ABF PP.H4, this study)
#   col 2 = Gobeil protein-QTL colocalization (PPH4)
#   right = derived interpretation (class + how the protein layer complements).
# Reads the frozen, verified join (build_external_pqtl_annotation.R). EXTERNAL,
# cited annotation — NOT re-derived, NOT a convergence-score input.
#
# Conventions: PDF, 6pt, all TEXT black, gene names italic, no title/subtitle
# (caption via message()); gradient fill light->magenta so in-cell numbers stay
# readable in black. Companion source CSV saved alongside.
# Output: figures/main/fig2_genetics/panels/FigS2O_pqtl_external_complementation.pdf
# Env: rnaseq
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/fig2_promoted_coloc_context.R"))

ann <- fread(file.path(BASE,
  "RNA-seq/results/multi_evidence/external_pqtl/pqtl_external_annotation.tsv"))
num <- function(x) suppressWarnings(as.numeric(x))

d <- ann[!is.na(num(gobeil_pqtl_coloc_pph4)) & gene != "DUSP23"]
d <- merge(d, load_fig2_promoted_context(BASE, d$gene), by = "gene", all.x = TRUE)
d[, expr := all_best]
d[, prot := num(gobeil_pqtl_coloc_pph4)]

# Derived interpretation from the promoted COLOC state plus the source-owned
# protein-QTL result. Missing or unsupported expression-QTL evidence is not
# described as a biological null.
d[, interp := fifelse(layer_relationship == "direction_discordant",
        "Expression-QTL support; eQTL/pQTL effect signs oppose",
     fifelse(all_support == "multi-signal COLOC",
        "Multi-signal expression-QTL support",
     fifelse(all_support == "single-signal COLOC only",
        "Single-signal expression-QTL support only",
     fifelse(all_support == "not evaluable",
        "Expression-QTL map not evaluable",
     fifelse(prot >= 0.80,
        "No PP.H4 > 0.5 support in evaluated expression-QTL map",
        "Protein-QTL below its 0.80 prioritization gate")))))]

# gene order top->bottom: flagship miss -> borderline -> sensitivity -> discordant -> reactive
ord <- c("MTARC1", "HSD17B13", "ERLIN1", "MTTP", "GCKR", "LGALS1")
d <- d[match(ord, gene)]
d[, gene := factor(gene, levels = rev(ord))]

fmt <- function(v) fifelse(v >= 0.995, "1.00", formatC(v, digits = 2, format = "fg"))
tiles <- rbindlist(list(
  d[, .(gene, xk = 1L, xlab = "Expression-QTL\n(this study)", value = expr, lab = fmt(expr))],
  d[, .(gene, xk = 2L, xlab = "Protein-QTL\n(Gobeil 2026)",  value = prot, lab = fmt(prot))]))

TXT_X <- 2.75
p <- ggplot() +
  geom_tile(data = tiles, aes(xk, gene, fill = value), colour = "white", linewidth = 0.6) +
  geom_text(data = tiles, aes(xk, gene, label = lab), size = GEOM_TEXT_6PT, colour = "black") +
  geom_text(data = d, aes(TXT_X, gene, label = interp), hjust = 0,
            size = GEOM_TEXT_6PT, colour = "black") +
  scale_fill_gradient(low = "#f3eef2", high = "#c0508d", limits = c(0, 1),
                      breaks = c(0, 0.5, 1), name = "coloc PP.H4") +
  scale_x_continuous(breaks = c(1, 2), labels = c("Expression-QTL\n(this study)",
                     "Protein-QTL\n(Gobeil 2026)"), limits = c(0.5, 7.0),
                     expand = expansion(mult = 0)) +
  labs(x = NULL, y = NULL) +
  theme_masld(base_size = 6) +
  theme(axis.text.y = element_text(face = "italic"),
        axis.line = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.key.size = unit(3, "mm"), legend.position = "bottom")

out_dir <- Sys.getenv("FIG2_SUPP_OUT_DIR",
  unset = file.path(BASE, "figures/main/fig2_genetics/panels"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_pdf <- file.path(out_dir, "FigS2O_pqtl_external_complementation.pdf")
fwrite(d[, .(gene, expr_qtl_coloc_ours = expr, protein_qtl_coloc_gobeil = prot,
             multi_signal_pph4 = all_susie, single_signal_pph4 = all_abf,
             expression_qtl_support = all_support,
             direct_multi_signal_pph4 = direct_susie,
             direct_single_signal_pph4 = direct_abf,
             direct_support_state = direct_support,
             enzyme_multi_signal_pph4 = enzyme_susie,
             enzyme_single_signal_pph4 = enzyme_abf,
             enzyme_support_state = enzyme_support, layer_relationship,
             interpretation = interp)],
       file.path(out_dir, "FigS2O_pqtl_external_complementation_source.csv"))
save_fig(p, out_pdf, width = 5.4, height = 2.3)
message("Wrote ", out_pdf)
message("CAPTION: Source-reported liver protein-QTL colocalization (Gobeil et al. 2026) ",
        "compared with the promoted expression-QTL COLOC portfolio. Expression-QTL support is ",
        "reported as multi-signal, single-signal only, evaluated without PP.H4 > 0.5 support, ",
        "or not evaluable. Cited external annotation; not a Resource evidence-class input.")
