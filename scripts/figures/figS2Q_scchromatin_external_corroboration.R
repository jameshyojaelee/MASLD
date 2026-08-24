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
source(file.path(BASE, "scripts/figures/fig2_promoted_coloc_context.R"))

src <- fread(file.path(BASE,
  "RNA-seq/results/multi_evidence/external_scchromatin/figS2Q_source.tsv"))
num <- function(v) suppressWarnings(as.numeric(v))
src <- merge(src, load_fig2_promoted_context(BASE, src$gene), by = "gene", all.x = TRUE)
deg <- load_fig2_current_deg_by_symbol(BASE, src$gene)
src <- merge(src, deg[, .(gene = gene_symbol, current_deg_ensembl = gene,
                          current_logFC = logFC, current_padj = padj,
                          current_is_deg)],
             by = "gene", all.x = TRUE)

# ---- rows to feature (grouped; ONECUT1/CEBPA/XBP1 stay in the data table only) --
ord <- c("RORA","PPP1R3B","KRT8","EFHD1","CUX2","RELB","CREB5","THRB","HNF4A",
         "TRPS1","FCGR2B","HLA-DQA1","FLACC1",
         "PNPLA3","TM6SF2","HSD17B13","MBOAT7")
G1 <- "Source comparison\n(regulatory layer)"
G2 <- "Cell-context\nannotation\n(their caQTL)"
G3 <- "Not discussed\nin source atlas"
grp <- c(rep(G1, 9), rep(G2, 4), rep(G3, 4))
names(grp) <- ord
d <- src[match(ord, gene)]
d[, group := factor(grp[gene], levels = c(G1, G2, G3))]

# ---- our genetic tiles: direct-MASLD vs enzyme-trait max PP.H4 ------------------
d[, gdir := direct_best]
d[, genz := enzyme_best]

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
d[, deg_txt := fifelse(is.na(current_is_deg), "bulk map not evaluable",
    fifelse(current_is_deg,
      sprintf("DEG %s (%.0e)", fifelse(num(current_logFC) > 0, "↑", "↓"),
              num(current_padj)), "no canonical DEG call"))]

# ---- release-derived interpretation retained in the source table --------------
d[, itext := sprintf("external: %s; direct/liver-fat: %s; enzyme: %s; bulk: %s",
                     elab, direct_support, enzyme_support, deg_txt)]

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
  geom_text(data = tiles, aes(xk, gene, label = lab), size = GEOM_TEXT_6PT, colour = "black") +
  geom_text(data = d, aes(TXT_EL, gene, label = elab), hjust = 0,
            size = GEOM_TEXT_6PT, colour = "black") +
  geom_text(data = d, aes(TXT_DEG, gene, label = deg_txt), hjust = 0,
            size = GEOM_TEXT_6PT, colour = "black") +
  scale_fill_gradient(low = "#f3eef2", high = "#c0508d", limits = c(0, 1),
                      breaks = c(0, 0.5, 1), name = "our coloc PP.H4") +
  scale_x_continuous(breaks = c(1, 2),
                     labels = c("direct /\nliver fat", "liver\nenzyme"),
                     limits = c(0.5, 7.0), expand = expansion(mult = 0),
                     sec.axis = sec_axis(~ ., breaks = c(1.5, TXT_EL + 0.35, TXT_DEG + 0.35),
                       labels = c("our coloc", "Elison sc-chromatin", "our DEG"))) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(x = "promoted expression-QTL COLOC", y = NULL) +
  theme_masld(base_size = 6) +
  theme(axis.text.y = element_text(face = "italic"),
        axis.line = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(),
        axis.title.x = element_text(size = 6, hjust = 0.08),
        strip.text.y.left = element_text(angle = 0, hjust = 0.5, size = 6),
        strip.placement = "outside", panel.spacing = unit(2, "mm"),
        legend.key.size = unit(3, "mm"), legend.position = "bottom")

out_dir <- Sys.getenv("FIG2_SUPP_OUT_DIR",
  unset = file.path(BASE, "figures/main/fig2_genetics/panels"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_pdf <- file.path(out_dir, "FigS2Q_scchromatin_external_corroboration.pdf")
fwrite(d[, .(gene, group, elison_layer = elab, our_genetic_direct_pph4 = gdir,
             our_genetic_enzyme_pph4 = genz, our_disease_state = deg_txt,
             current_bulk_ensembl = current_deg_ensembl,
             current_bulk_logFC = current_logFC,
             current_bulk_padj = current_padj,
             current_bulk_is_deg = current_is_deg,
             elison_finding_type, elison_celltype,
             direct_multi_signal_pph4 = direct_susie,
             direct_single_signal_pph4 = direct_abf,
             direct_support_state = direct_support,
             enzyme_multi_signal_pph4 = enzyme_susie,
             enzyme_single_signal_pph4 = enzyme_abf,
             enzyme_support_state = enzyme_support,
             interpretation = itext)],
       file.path(out_dir, "FigS2Q_scchromatin_external_corroboration_source.csv"))
save_fig(p, out_pdf, width = 6.3, height = 3.4)
message("Wrote ", out_pdf)
message("CAPTION: Source-reported single-cell chromatin annotations from Elison/Gaulton et al. ",
        "are cross-referenced to the promoted expression-QTL COLOC portfolio and the current ",
        "canonical bulk gate. Direct/liver-fat and liver-enzyme posteriors remain separate. ",
        "The source table records method availability and support states; absent support is not ",
        "treated as proof against a regulatory mechanism. External annotation; not a Resource ",
        "evidence-class input.")
