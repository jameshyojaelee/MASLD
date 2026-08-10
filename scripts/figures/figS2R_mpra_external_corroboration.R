#!/usr/bin/env Rscript
# figS2R_mpra_external_corroboration.R
# ─────────────────────────────────────────────────────────────────────────────
# Fig 2 supplement (FigS2R): external EXPERIMENTAL variant-function corroboration
# of our genetic + disease-state maps, from the MASLD MPRA + sc-CRISPRi study of
# Zhu / Hu et al. 2026 (Nat Genet; DOI 10.1038/s41588-026-02617-8). Compact
# tile-table:
#   tile 1 = our genetic map, DIRECT-MASLD colocalization (max SuSiE/ABF PP.H4)
#   tile 2 = our genetic map, ENZYME-trait colocalization (max SuSiE/ABF PP.H4)
#            -- side by side to keep trait-anchoring transparent
#   text   = Zhu MPRA layer (cell model + allelic direction) | our DEG | verdict
# Rows grouped: featured targets | corroborated DAVs (our genetic/disease map) |
# expression-miss (Zhu MPRA-functional where our expression-genetics is null).
# EXTERNAL, cited annotation -- NOT re-derived, NOT a convergence-score input.
#
# Conventions: PDF, 6pt, all TEXT black, gene names italic, no title/subtitle
# (caption via message()); light->magenta fill so in-cell numbers stay readable.
# Reads the frozen join (build_external_mpra_annotation.R).
# Out: figures/main/fig2_genetics/panels/FigS2R_mpra_external_corroboration.pdf
# Env: rnaseq
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))

src <- fread(file.path(BASE,
  "RNA-seq/results/multi_evidence/external_mpra/figS2R_source.tsv"))
num <- function(v) suppressWarnings(as.numeric(v))

# ---- our genetic tiles: direct-MASLD vs enzyme-trait max PP.H4 ------------------
src[, gdir := pmax(fcoalesce(num(our_susie_direct),0), fcoalesce(num(our_abf_direct),0))]
src[, genz := pmax(fcoalesce(num(our_susie_enzyme),0), fcoalesce(num(our_abf_enzyme),0))]

# ---- row groups (featured always top; then corroborated; then gaps) ------------
G1 <- "Featured\ntargets"
G2 <- "Corroborated\n(our genetic /\ndisease map)"
G3 <- "Expression-miss\n(their MPRA,\nour maps null)"
src[, group := fifelse(zhu_finding_type == "featured_target", G1,
                fifelse(relationship == "expression_miss", G3, G2))]
src[, group := factor(group, levels = c(G1, G2, G3))]

# order within group: featured by curated order; others by direct then enzyme coloc
feat_ord <- c("SLC22A3","LPL","ANGPTL3","APOA5")
src[, ford := fifelse(group == G1, match(gene, feat_ord), NA_integer_)]
setorder(src, group, ford, -gdir, -genz, gene)
ord <- src$gene
src[, gene := factor(gene, levels = rev(ord))]

# ---- Zhu MPRA layer label (cell model + allelic direction) ---------------------
arr <- function(d) fifelse(d == "activity_up", "↑",
                    fifelse(d == "activity_down", "↓", ""))
src[, zlab := fifelse(zhu_finding_type == "featured_target", "featured DAV target",
               sprintf("%s MPRA %s (log2FC %.1f)", zhu_celltype, arr(zhu_direction),
                       num(zhu_mpra_log2fc)))]

# ---- our disease-state DEG text ------------------------------------------------
src[, deg_txt := fifelse(our_is_deg == TRUE,
    sprintf("DEG %s (%.0e)", fifelse(num(our_bulk_logFC) > 0, "↑", "↓"),
            num(our_bulk_treat_fdr)), "n.s.")]

# ---- short verdict -------------------------------------------------------------
vlab <- c(corroborated_genetic_direct = "direct-MASLD coloc",
          corroborated_genetic_enzyme = "enzyme-trait coloc",
          corroborated_disease_state  = "disease-state DEG",
          expression_miss             = "gap: MPRA adds signal")
src[, vtext := vlab[relationship]]

# ---- layout --------------------------------------------------------------------
fmt <- function(v) fifelse(v >= 0.995, "1.00",
                    fifelse(v < 0.005, "–", formatC(v, digits = 2, format = "fg")))
tiles <- rbindlist(list(
  src[, .(gene, group, xk = 1L, value = gdir, lab = fmt(gdir))],
  src[, .(gene, group, xk = 2L, value = genz, lab = fmt(genz))]))

TXT_Z <- 2.85; TXT_DEG <- 5.35; TXT_V <- 7.15
p <- ggplot() +
  geom_tile(data = tiles, aes(xk, gene, fill = value), colour = "white", linewidth = 0.5) +
  geom_text(data = tiles, aes(xk, gene, label = lab), size = 1.9, colour = "black") +
  geom_text(data = src, aes(TXT_Z,   gene, label = zlab),    hjust = 0, size = 1.9, colour = "black") +
  geom_text(data = src, aes(TXT_DEG, gene, label = deg_txt), hjust = 0, size = 1.9, colour = "black") +
  geom_text(data = src, aes(TXT_V,   gene, label = vtext),   hjust = 0, size = 1.9, colour = "black") +
  scale_fill_gradient(low = "#f3eef2", high = "#c0508d", limits = c(0, 1),
                      breaks = c(0, 0.5, 1), name = "our coloc PP.H4") +
  scale_x_continuous(breaks = c(1, 2),
                     labels = c("direct\nMASLD", "enzyme\ntrait"),
                     limits = c(0.5, 9.2), expand = expansion(mult = 0),
                     sec.axis = sec_axis(~ ., breaks = c(1.5, TXT_Z + 0.3, TXT_DEG + 0.3, TXT_V + 0.3),
                       labels = c("our coloc", "Zhu MPRA", "our DEG", "verdict"))) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(x = "our genetic map (colocalization)", y = NULL) +
  theme_masld(base_size = 6) +
  theme(axis.text.y = element_text(face = "italic"),
        axis.line = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(),
        axis.title.x = element_text(size = 5.5, hjust = 0.06),
        strip.text.y.left = element_text(angle = 0, hjust = 0.5, size = 5),
        strip.placement = "outside", panel.spacing = unit(2, "mm"),
        legend.key.size = unit(3, "mm"), legend.position = "bottom")

out_dir <- file.path(BASE, "figures/main/fig2_genetics/panels")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_pdf <- file.path(out_dir, "FigS2R_mpra_external_corroboration.pdf")
fwrite(src[, .(gene, group, zhu_mpra = zlab, our_genetic_direct_pph4 = gdir,
               our_genetic_enzyme_pph4 = genz, our_disease_state = deg_txt,
               relationship, verdict = vtext, caveat)],
       file.path(out_dir, "FigS2R_mpra_external_corroboration_source.csv"))
save_fig(p, out_pdf, width = 5.5, height = 0.13 * nrow(src) + 1.4)
message("Wrote ", out_pdf)
message("CAPTION: External experimental variant-function corroboration (Zhu, Hu et al. 2026, ",
        "Nat Genet; MASLD MPRA in HepG2/LX-2 + sc-CRISPRi). Tiles = our colocalization PP.H4 ",
        "anchored on direct-MASLD vs liver-enzyme traits (side by side to show trait anchoring); ",
        "text = Zhu's MPRA differential-activity call (cell model + allelic direction), our ",
        "canonical DEG, and verdict. Featured targets: LPL corroborates our disease-state map ",
        "(DEG logFC +2.8), SLC22A3 our enzyme-anchored genetic map; APOA5 and ANGPTL3 are ",
        "expression-miss gaps where the MPRA flags cis-regulatory function our expression-genetics ",
        "is null on. MPRA/CRISPRi are cell-line reporter assays (not primary-tissue), and the gene ",
        "label on a DAV row is our finemapping/coloc locus assignment (the DAV and log2FC are ",
        "Zhu's). Cited external annotation, not a convergence-score input.")
