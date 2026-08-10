#!/usr/bin/env Rscript
# ==============================================================================
# Fig 4d — CROSS-LINEAGE open-chromatin context for prioritized targets.
# KEY MESSAGE: Fine-mapped liver-trait variants intersect open chromatin across several liver lineages but provide regulatory context rather than functional validation.
# Reframe of the retired hepatocyte-only panel (archived under
# _supp/snatac_accessibility_hepatocyte_only_superseded.pdf). The claim is SCOPE,
# not a single lineage: fine-mapped GWAS credible-set variants of prioritized
# MASLD targets sit in OPEN snATAC peaks (GSE244832) distributed ACROSS cell
# lineages — hepatocytes are one of several, essentially tied with fibroblasts —
# yet only a small minority ALSO colocalize with an eQTL. Broad regulatory
# potential, narrow colocalization. This is STATIC accessibility (variant in an
# open peak), NOT disease-state differential accessibility (that GSE244832
# contrast is a 5-control, caQTL-anti-concordant liability → excluded).
#
#   Left  (bars)     : per lineage, N prioritized targets with a variant in an
#                      open peak (gray) with the also-colocalizing subset (red).
#   Right (beeswarm) : the full prioritized set — x = fine-map PIP; fill = COLOC
#                      PP.H4 (PolyFun); ring = variant disrupts a TF motif
#                      (motifbreakR, putative). Breadth lives in the left bars.
#
# Output: figures/main/fig4_validation/panels/fig4d_snatac_accessibility.pdf
# Env: rnaseq
# ==============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrepel); library(patchwork) })
set.seed(42)
BASE <- Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
# Non-canonical / artifact-looking symbols dropped in the sibling hepatocyte panel;
# gene-level, so they apply equally to the broader cross-lineage set.
DROP <- c("SPECC1L-ADORA2A","RPS24","EXOC3L4","LRRC75B","SLC50A1","CKM","DYNLRB2","MAU2",
          "GATAD2A","RFXANK","SPPL3","FLAD1","ZNF827","SPATA31H1","SIPA1L2","CWF19L1","FCHO2")
# Consolidate labels to include only key, established, or prioritized MASLD targets
# to avoid crowded overlap (hecticness) in the panel's Fig4 layout slot (3.13x1.88in).
LAB <- c("CYP26A1", "GCKR", "GPAM", "MTTP", "P2RX7")
# Lineage groups: the immune lymphoid compartment (T + resident NK + plasma) is
# consolidated into a single "Lymphoid" bar; counts are UNIONS (a gene open in
# more than one member lineage is still counted once).
GRP <- list(
  Hepatocytes    = "Hepatocytes",
  Fibroblasts    = "Fibroblasts",
  Macrophages    = "Macrophages",
  Cholangiocytes = "Cholangiocytes",
  Lymphoid       = c("T_cells", "Resident_NK", "Plasma_cells"),
  Endothelial    = "Endothelial_cells")

# ---- load + filter (identical universe / biotype / DROP gates to the sibling panel) ----
d <- fread(file.path(BASE, "GWAS/finemapping/results/gwas_atac/gene_level_gwas_atac.csv"))
uni <- trimws(readLines(file.path(BASE, "Analysis/Spatial/results/universe_validation/prioritized_universe_FINAL.txt")))
uni <- uni[uni != ""]
d <- d[gwas_variant_in_peak == TRUE & assigned_gene %in% uni]
atl <- fread(file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"))
gcol <- if ("human_symbol" %in% names(atl)) "human_symbol" else "gene"
num <- function(x) suppressWarnings(as.numeric(x))
atl[, pp4 := pmax(num(coloc_best_susie_pp4_polyfun), num(coloc_best_pp4_polyfun), na.rm = TRUE)]
am <- unique(atl[, .(gene = get(gcol), pp4, gene_biotype)], by = "gene")
d <- merge(d[, .(gene = assigned_gene, cts = gwas_variant_cell_types,
                 pip = as.numeric(gwas_max_pip_in_peak),
                 motif = !is.na(gwas_motif_disrupted) & gwas_motif_disrupted != "")],
           am, by = "gene", all.x = TRUE)
d[is.na(pp4), pp4 := 0]
# funnel counts for the caption (nesting within the shared 9,882-gene universe):
n_inpeak_any <- length(unique(d$gene))                                   # variant-in-peak ∩ universe, any biotype
n_inpeak_pc  <- length(unique(d[gene_biotype == "protein_coding"]$gene)) # + protein-coding
d <- d[gene_biotype == "protein_coding" & !(gene %in% DROP)]
d[, coloc := pp4 > 0.5]

# ---- per-lineage set sizes (a gene can be open in several lineages) ----
res <- rbindlist(lapply(names(GRP), function(g){
  pat <- paste0("(", paste(GRP[[g]], collapse = "|"), ")")   # union match across group members
  sub <- d[grepl(pat, cts)]
  data.table(lineage = g, n_genes = nrow(sub), n_coloc = sum(sub$coloc))
}))
res[, disp := gsub("_", " ", lineage)]
res[, disp := factor(disp, levels = disp[order(n_genes)])]
cat(sprintf("[fig4d-snatac] %d prioritized protein-coding targets, variant in an open peak; colocalizing (PP.H4>0.5)=%d; motif=%d\n",
            nrow(d), sum(d$coloc), sum(d$motif)))
print(res[order(-n_genes)])

# ---- (i) per-lineage bullet bars: breadth (gray) + narrow colocalization (red) ----
pa <- ggplot(res, aes(y = disp)) +
  geom_col(aes(x = n_genes, fill = "in open peak"), width = 0.68) +
  geom_col(aes(x = n_coloc, fill = "+ eQTL-colocalizes"), width = 0.68) +
  geom_text(aes(x = n_genes, label = n_genes), hjust = -0.3, size = 6/.pt, colour = house_ink) +
  scale_fill_manual(values = c("in open peak" = "#BDBDBD", "+ eQTL-colocalizes" = "#B2182B"),
                    breaks = c("in open peak", "+ eQTL-colocalizes"), name = NULL,
                    guide = guide_legend(ncol = 1)) +
  scale_x_continuous(limits = c(0, max(res$n_genes) * 1.16), breaks = c(0, 30, 60), expand = c(0, 0)) +
  labs(x = "targets, open-peak variant", y = NULL) +
  theme_masld() +
  theme(legend.position = "bottom", legend.direction = "vertical",
        legend.key.size = unit(0.3, "lines"),
        legend.text = element_text(size = 6), legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(2, "pt"), panel.grid = element_blank(),
        axis.ticks.y = element_blank(), axis.line.y = element_blank())

# ---- (ii) beeswarm of the full set: PIP x COLOC-fill x motif-ring x breadth-size ----
setorder(d, pip)
d[, yj := ((seq_len(.N) %% 9) - 4) / 3.4]
d[, yj := yj + runif(.N, -0.12, 0.12)]
lab <- d[gene %in% LAB]
pb <- ggplot(d, aes(x = pip, y = yj)) +
  geom_point(aes(fill = pp4, colour = motif), shape = 21, size = 1.1, stroke = 0.35) +
  ggrepel::geom_text_repel(data = lab, aes(x = pip, y = yj, label = gene), size = 6/.pt, fontface = "italic",
                           nudge_x = ifelse(lab$gene == "GCKR", -0.15, 0),
                           colour = house_ink, max.overlaps = 40, min.segment.length = 0, segment.size = 0.2,
                           segment.colour = "grey70", box.padding = 0.75, point.padding = 0.35, force = 6,
                           max.iter = 10000, max.time = 2, seed = 7) +
  scale_fill_gradient(low = "#ECECEC", high = "#B2182B", limits = c(0, 1), breaks = c(0, 0.5, 1), name = "COLOC PP.H4",
                      guide = guide_colourbar(direction = "horizontal", title.position = "top", title.hjust = 0.5,
                                              barwidth = unit(1.6, "cm"), barheight = unit(0.16, "cm"), order = 1)) +
  scale_colour_manual(values = c(`TRUE` = house_ink, `FALSE` = "grey65"),
                      breaks = c(TRUE, FALSE), labels = c("disrupts TF motif", "no motif hit"), name = NULL,
                      guide = guide_legend(order = 2, ncol = 1, override.aes = list(size = 2))) +
  scale_x_continuous(limits = c(-0.02, 1.04), breaks = c(0, 0.5, 1.0), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-1.7, 1.7), breaks = NULL, expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(x = "fine-map PIP (open-peak variant)", y = NULL) +
  theme_masld() +
  theme(legend.position = "bottom", legend.box = "horizontal", legend.key.size = unit(0.3, "lines"),
        legend.text = element_text(size = 6), legend.title = element_text(size = 6),
        legend.margin = margin(0, 0, 0, 0), legend.spacing.x = unit(4, "pt"),
        panel.grid = element_blank(), axis.ticks.y = element_blank())

p <- pa + pb + plot_layout(widths = c(0.8, 1.9))
out <- file.path(FIG4_DIR, "panels", "fig4d_snatac_accessibility.pdf")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
ggsave(out, p, width = 3.13, height = 1.88, device = grDevices::cairo_pdf)
cat("[fig4d-snatac] saved:", out, "\n")

message(sprintf(paste0("CAPTION (Fig 4d): Cross-lineage open-chromatin context for prioritized targets. ",
  "From the current %s-gene prioritized universe, ",
  "%d targets carry a fine-mapped credible-set variant inside an open snATAC peak (any biotype); restricting to ",
  "protein-coding (%d) and removing %d artifact-looking symbols yields the set shown. ",
  "Of %d canonical protein-coding prioritized MASLD targets whose fine-mapped GWAS credible-set variant lies inside an ",
  "accessible snATAC peak (GSE244832), the open-peak lineage is broadly distributed — Hepatocytes %d, Fibroblasts %d, ",
  "Macrophages %d, Cholangiocytes %d, Lymphoid (T/NK/plasma) %d, Endothelial %d — NOT hepatocyte-exclusive; ",
  "hepatocytes are essentially tied with fibroblasts. Left: per-lineage counts (gray) with the also-eQTL-colocalizing subset ",
  "(red, PP.H4>0.5, PolyFun); only %d of %d targets colocalize in any lineage. Right: the full set — x = fine-map PIP of the ",
  "credible-set variant; fill = COLOC PP.H4; black ring = the variant also disrupts a TF motif (motifbreakR, PUTATIVE, not ",
  "validated). IMPORTANT: this is STATIC accessibility (a variant sits ",
  "in an open peak), NOT disease-state differential accessibility — the GSE244832 disease-DA contrast (5 controls, caQTL ",
  "anti-concordant, p=1.6e-4) is an established liability and is EXCLUDED. The gene-level display does not require PIP, ",
  "PP.H4, and motif status to arise from the same variant/trait locus. Peaks are accessible-IN-a-lineage, not lineage-EXCLUSIVE ",
  "(a gene can be open across several lineages). The panel is regulatory context rather than functional validation; its ",
  "statistics are not re-derived here."),
  format(length(uni), big.mark=","), n_inpeak_any, n_inpeak_pc, n_inpeak_pc - nrow(d),
  nrow(d), res[lineage=="Hepatocytes", n_genes], res[lineage=="Fibroblasts", n_genes],
  res[lineage=="Macrophages", n_genes], res[lineage=="Cholangiocytes", n_genes],
  res[lineage=="Lymphoid", n_genes], res[lineage=="Endothelial", n_genes],
  sum(d$coloc), nrow(d)))
