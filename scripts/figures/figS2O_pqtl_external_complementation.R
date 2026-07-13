#!/usr/bin/env Rscript
# figS2O_pqtl_external_complementation.R
# ─────────────────────────────────────────────────────────────────────────────
# Fig 2 supplement (FigS2O): external liver protein-QTL complementation of our
# expression-QTL genetic map. Compact heatmap-table for the liver-proteogenomic
# causal effectors of Gobeil et al. 2026 (medRxiv 10.64898/2026.06.04.26354903):
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

ann <- fread(file.path(BASE,
  "RNA-seq/results/multi_evidence/external_pqtl/pqtl_external_annotation.tsv"))
num <- function(x) suppressWarnings(as.numeric(x))

d <- ann[!is.na(num(gobeil_pqtl_coloc_pph4)) & gene != "DUSP23"]
d[, expr := pmax(fcoalesce(num(our_max_susie_pp4), 0), fcoalesce(num(our_max_abf_pp4), 0))]
d[, prot := num(gobeil_pqtl_coloc_pph4)]

# derived interpretation (class + how protein-QTL complements) — from data fields
short_class <- function(c) fifelse(c == "genetic_only", "Genetic-only",
                          fifelse(c == "disease_state_only", "Disease-state-only", "Neither"))
d[, interp := fifelse(layer_relationship == "direction_discordant",
        "Genetic-only · colocalizes on expression; eQTL/pQTL effect signs oppose",
     fifelse(our_primary_class == "disease_state_only",
        "Disease-state-only · reactive in both maps (their analysis agrees)",
     fifelse(gobeil_pqtl_hit & our_genetic_sensitivity == TRUE,
        "Neither (ABF-suggestive) · protein-QTL confirms the genetic signal",
     fifelse(gobeil_pqtl_hit & our_genetic_sensitivity == FALSE,
        "Neither · expression-map miss; protein-QTL recovers it",
        "Neither · expression-map miss; protein-QTL borderline (<0.80)"))))]

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
  geom_text(data = tiles, aes(xk, gene, label = lab), size = 2.1, colour = "black") +
  geom_text(data = d, aes(TXT_X, gene, label = interp), hjust = 0, size = 2.0, colour = "black") +
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

out_dir <- file.path(BASE, "figures/main/fig2_genetics/panels")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_pdf <- file.path(out_dir, "FigS2O_pqtl_external_complementation.pdf")
fwrite(d[, .(gene, expr_qtl_coloc_ours = expr, protein_qtl_coloc_gobeil = prot,
             our_primary_class, layer_relationship, interpretation = interp)],
       file.path(out_dir, "FigS2O_pqtl_external_complementation_source.csv"))
save_fig(p, out_pdf, width = 5.4, height = 2.3)
message("Wrote ", out_pdf)
message("CAPTION: External protein-QTL complementation (Gobeil et al. 2026, liver mass-spec pQTL) ",
        "of our expression-QTL map. Of their five colocalized effectors we already flag four ",
        "(GCKR primary; ERLIN1/MTTP sensitivity; LGALS1 reactive) and miss one (MTARC1); GCKR is ",
        "eQTL/pQTL direction-discordant. Cited external annotation, not a convergence-score input.")
