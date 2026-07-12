#!/usr/bin/env Rscript
# figS2Q_scchromatin_external_corroboration.R
# ─────────────────────────────────────────────────────────────────────────────
# Fig 2 supplement (FigS2Q): external single-cell-chromatin corroboration of our
# genetic + disease-state maps, and the complementary-layer boundary, from the
# liver single-cell multiomics atlas of Elison/Gaulton et al. 2025
# (medRxiv 10.1101/2025.05.09.25327043). Compact tile-table:
#   tile 1 = our genetic map, DIRECT-MASLD colocalization (max SuSiE/ABF PP.H4)
#   tile 2 = our genetic map, ENZYME-trait colocalization (max SuSiE/ABF PP.H4)
#            -- the two side-by-side make the trait-anchoring transparent
#   text   = Elison's single-cell chromatin layer | our disease-state DEG | interpretation
# Rows grouped: corroboration (their GRN/coloc vs our maps) | cell-of-action
# gap-fill (their per-cell-type caQTL) | mutual boundary (lipid canon absent both).
# EXTERNAL, cited annotation -- NOT re-derived, NOT a convergence-score input.
#
# Conventions: PDF, 6pt, all TEXT black, gene names italic, no title/subtitle
# (caption via message()); light->magenta fill so in-cell numbers stay readable.
# Reads the frozen join (build_external_scchromatin_annotation.R).
# Out: figures/main/fig2_genetics/panels/FigS2Q_scchromatin_external_corroboration.pdf
# Env: rnaseq
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))

src <- fread(file.path(BASE,
  "RNA-seq/results/multi_evidence/external_scchromatin/figS2Q_source.tsv"))
num <- function(v) suppressWarnings(as.numeric(v))

# ---- rows to feature (grouped; ONECUT1/CEBPA/XBP1 stay in the data table only) --
ord <- c("RORA","PPP1R3B","KRT8","EFHD1","CUX2","RELB","CREB5","THRB","HNF4A",
         "TRPS1","FCGR2B","HLA-DQA1","FLACC1",
         "PNPLA3","TM6SF2","HSD17B13","MBOAT7")
G1 <- "Corroboration\n(their regulatory\nlayer vs our maps)"
G2 <- "Cell-of-action\ngap-fill\n(their caQTL)"
G3 <- "Mutual boundary\n(lipid canon)"
grp <- c(rep(G1, 9), rep(G2, 4), rep(G3, 4))
names(grp) <- ord
d <- src[match(ord, gene)]
d[, group := factor(grp[gene], levels = c(G1, G2, G3))]

# ---- our genetic tiles: direct-MASLD vs enzyme-trait max PP.H4 ------------------
d[, gdir := pmax(fcoalesce(num(our_susie_direct),0), fcoalesce(num(our_abf_direct),0))]
d[, genz := pmax(fcoalesce(num(our_susie_enzyme),0), fcoalesce(num(our_abf_enzyme),0))]

# ---- Elison single-cell chromatin layer label ---------------------------------
elison_lab <- c(RORA="MASL hep GRN ↑", HNF4A="MASL hep GRN ↑", CUX2="MASL hep GRN ↑",
  CREB5="Fib+ hep GRN ↑", RELB="Fib+ hep GRN ↑", THRB="hep GRN ↓ (drug axis)",
  EFHD1="hep eQTL-coloc", KRT8="hep loop (no eQTL)", PPP1R3B="hep caQTL + loop >200kb",
  TRPS1="myeloid caQTL", FCGR2B="endothelial caQTL", "HLA-DQA1"="myeloid caQTL",
  FLACC1="HSC caQTL",
  PNPLA3="absent (0 mentions)", TM6SF2="absent (0 mentions)",
  HSD17B13="absent (0 mentions)", MBOAT7="absent (0 mentions)")
d[, elab := elison_lab[gene]]

# ---- our disease-state DEG text ------------------------------------------------
d[, deg_txt := fifelse(our_is_deg == TRUE,
    sprintf("DEG %s (%.0e)", fifelse(num(our_bulk_logFC) > 0, "↑", "↓"),
            num(our_bulk_treat_fdr)),
    fifelse(gene %in% c("CREB5"), "near-DEG (0.057)", "n.s."))]

# ---- interpretation (short; honest) -------------------------------------------
interp <- c(
  RORA   = "genetic (direct 0.98) + DEG + their GRN — co-nomination (axes differ)",
  PPP1R3B= "genetic (direct-MASLD ABF 0.97); their causal variant + >200 kb loop add cell type",
  KRT8   = "we supply expression-COLOC they lacked; enzyme-trait anchored",
  EFHD1  = "co-localized locus; direction differs (their eQTL ↑ vs our DE ↓)",
  CUX2   = "disease-state DEG + suggestive coloc; their MASL GRN ↑ / our pooled DE ↓",
  RELB   = "disease-state DEG ↑ matches their Fib+ GRN ↑ (same axis)",
  CREB5  = "our Progressor-subtype marker + Cas13 target = their MASH GRN driver",
  THRB   = "drug-axis: their GRN ↓ + our regulon ↓; GGT-anchored, direct null",
  HNF4A  = "TF-activity agreement; our expression-genetic null (stable master TF)",
  TRPS1  = "their myeloid caQTL assigns the cell of action for our genetic hit",
  FCGR2B = "their endothelial caQTL assigns the cell of action for our genetic hit",
  "HLA-DQA1"="their myeloid caQTL resolves the cell layer",
  FLACC1 = "their HSC caQTL detects signal our expression-COLOC is null on",
  PNPLA3 = "held by monogenic I148M — blind to caQTL AND expression-COLOC",
  TM6SF2 = "held by coding E167K — direct-MASLD coloc null",
  HSD17B13="held by splice rs72613567 — expression-COLOC null",
  MBOAT7 = "rescued by TWAS (z = -3.69) — coloc null both maps")
d[, itext := interp[gene]]

# ---- layout: gene on y (reverse group order top->bottom) -----------------------
d[, gene := factor(gene, levels = rev(ord))]
fmt <- function(v) fifelse(v >= 0.995, "1.00",
                    fifelse(v < 0.005, "–", formatC(v, digits = 2, format = "fg")))
tiles <- rbindlist(list(
  d[, .(gene, group, xk = 1L, value = gdir, lab = fmt(gdir))],
  d[, .(gene, group, xk = 2L, value = genz, lab = fmt(genz))]))

TXT_EL <- 2.85; TXT_DEG <- 5.6
p <- ggplot() +
  geom_tile(data = tiles, aes(xk, gene, fill = value), colour = "white", linewidth = 0.6) +
  geom_text(data = tiles, aes(xk, gene, label = lab), size = 2.0, colour = "black") +
  geom_text(data = d, aes(TXT_EL,  gene, label = elab),   hjust = 0, size = 2.0, colour = "black") +
  geom_text(data = d, aes(TXT_DEG, gene, label = deg_txt),hjust = 0, size = 2.0, colour = "black") +
  scale_fill_gradient(low = "#f3eef2", high = "#c0508d", limits = c(0, 1),
                      breaks = c(0, 0.5, 1), name = "our coloc PP.H4") +
  scale_x_continuous(breaks = c(1, 2),
                     labels = c("direct\nMASLD", "enzyme\ntrait"),
                     limits = c(0.5, 7.0), expand = expansion(mult = 0),
                     sec.axis = sec_axis(~ ., breaks = c(1.5, TXT_EL + 0.35, TXT_DEG + 0.35),
                       labels = c("our coloc", "Elison sc-chromatin", "our DEG"))) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(x = "our genetic map (colocalization)", y = NULL) +
  theme_masld(base_size = 6) +
  theme(axis.text.y = element_text(face = "italic"),
        axis.line = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(),
        axis.title.x = element_text(size = 5.5, hjust = 0.08),
        strip.text.y.left = element_text(angle = 0, hjust = 0.5, size = 5),
        strip.placement = "outside", panel.spacing = unit(2, "mm"),
        legend.key.size = unit(3, "mm"), legend.position = "bottom")

out_dir <- file.path(BASE, "figures/main/fig2_genetics/panels")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_pdf <- file.path(out_dir, "FigS2Q_scchromatin_external_corroboration.pdf")
fwrite(d[, .(gene, group, elison_layer = elab, our_genetic_direct_pph4 = gdir,
             our_genetic_enzyme_pph4 = genz, our_disease_state = deg_txt,
             relationship, interpretation = itext, caveat)],
       file.path(out_dir, "FigS2Q_scchromatin_external_corroboration_source.csv"))
save_fig(p, out_pdf, width = 6.3, height = 3.4)
message("Wrote ", out_pdf)
message("CAPTION: External single-cell-chromatin corroboration (Elison, Gaulton et al. 2025; ",
        "liver 5-modality multiomics, 86 donors) of our genetic and disease-state maps. Tiles = our ",
        "colocalization PP.H4 anchored on direct-MASLD vs enzyme traits (side by side to show trait ",
        "anchoring); text = Elison's cell-type regulatory layer, our canonical DEG, and interpretation. ",
        "RORA/PPP1R3B corroborate our genetic map; RORA/CUX2/RELB/EFHD1/KRT8 our disease-state map; ",
        "per-cell-type caQTL resolves the cell-of-action for our genetic-only hits (TRPS1/FCGR2B/",
        "HLA-DQA1/FLACC1); the lipid canon (PNPLA3/TM6SF2/HSD17B13/MBOAT7) is absent from their atlas ",
        "AND null on our expression-COLOC, held by our TWAS/monogenic layers. Direction/axis and ",
        "trait-anchor caveats per gene in the source CSV. Cited external annotation, not a ",
        "convergence-score input; their genome-wide tables remain embargoed (Zenodo 15298484).")
