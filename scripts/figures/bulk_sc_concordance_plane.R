#!/usr/bin/env Rscript
# bulk_sc_concordance_plane.R  — CANDIDATE panel (evidence-plane grammar)
# Bulk RNA-seq logFC vs hepatocyte scRNA pseudobulk logFC, per gene. The honest
# per-gene view of the bulk<->scRNA concordance that celltype_concordance.R shows
# only as an aggregate Spearman rho (hep ~0.31). Mirrors coloc_deg_bridge.R.
#
# x = bulk logFC (canonical_deg_results.csv), y = hepatocyte pseudobulk logFC
# (MASLD_vs_Healthy). Gray cloud = all overlapping genes; y=x = perfect
# concordance; anchors = the Fig 2/3/4 cast, marks coloured concordant/discordant,
# TEXT black. Original celltype_concordance.R panel is UNTOUCHED.
#
# Output (PDF only): FIG4_SC_DIR/panels/bulk_sc_concordance_hep.pdf (exploratory)

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggrepel); library(ggrastr)
})
set.seed(42)
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))
PANEL_DIR <- file.path(FIG4_SC_DIR, "panels", "supplementary"); DATA_DIR <- file.path(PANEL_DIR, "data")
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
lab_size <- 6 / ggplot2::.pt

# ── Data: bulk logFC (symbol + ENSG) × hepatocyte pseudobulk logFC ────────────
bulk <- fread(file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv"),
  select = c("gene", "symbol", "logFC"))
setnames(bulk, "logFC", "bulk_lfc")
hep <- fread(file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2/pseudobulk_de/Hepatocytes_de.csv"),
  select = c("gene", "logFC", "contrast"))
stopifnot(all(hep$contrast == "MASLD_vs_Healthy"))
setnames(hep, "logFC", "sc_lfc")

# sc `gene` is a mix of HGNC symbols + ENSG — dual-key join (as celltype_concordance.R)
hep[, is_ensg := grepl("^ENSG", gene)]
m_sym  <- merge(hep[is_ensg == FALSE], bulk[, .(symbol, bulk_lfc)],
                by.x = "gene", by.y = "symbol")
m_ensg <- merge(hep[is_ensg == TRUE],  bulk[, .(gene, bulk_lfc)], by = "gene")
dt <- rbind(m_sym[, .(gene, bulk_lfc, sc_lfc)], m_ensg[, .(gene, bulk_lfc, sc_lfc)])
dt <- unique(dt[!is.na(bulk_lfc) & !is.na(sc_lfc)])

rho <- cor(dt$bulk_lfc, dt$sc_lfc, method = "spearman")
cat(sprintf("[bulk<->sc hep] n=%d genes, Spearman rho=%.3f\n", nrow(dt), rho))

# ── Anchors: Fig 2/3/4 cast; colour by concordance (same sign in bulk & sc) ───
anchors <- c("HKDC1","THRB","RORA","CYP3A4","AKR1B10","SERPINE1","FASN","SCD","FGF21")
dt[, hl := fifelse(gene %in% anchors,
                   fifelse(sign(bulk_lfc) == sign(sc_lfc), "concordant", "discordant"),
                   NA_character_)]
lab <- dt[gene %in% anchors]
COL_CONC <- "#00897B"   # teal — bulk & sc agree
COL_DISC <- "#C9265E"   # magenta — disagree
COL_BG   <- masld_colors$control
lim <- max(abs(c(dt$bulk_lfc, dt$sc_lfc)), na.rm = TRUE)

p <- ggplot() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.3, colour = "gray55") +
  geom_hline(yintercept = 0, linewidth = 0.2, colour = "gray85") +
  geom_vline(xintercept = 0, linewidth = 0.2, colour = "gray85") +
  rasterize(geom_point(data = dt, aes(bulk_lfc, sc_lfc),
                       colour = COL_BG, size = 0.35, alpha = 0.25, shape = 16), dpi = 600) +
  geom_point(data = dt[!is.na(hl)], aes(bulk_lfc, sc_lfc, colour = hl), size = 1.7, shape = 16) +
  geom_point(data = dt[!is.na(hl)], aes(bulk_lfc, sc_lfc),
             colour = "white", size = 1.7, shape = 1, stroke = 0.3) +
  geom_text_repel(data = lab, aes(bulk_lfc, sc_lfc, label = gene),
    colour = "black", size = lab_size, fontface = "italic",
    box.padding = 0.45, point.padding = 0.25, segment.size = 0.2, segment.color = "gray60",
    min.segment.length = 0, max.overlaps = Inf, seed = 42, force = 7,
    bg.color = "white", bg.r = 0.12) +
  scale_colour_manual(values = c(concordant = COL_CONC, discordant = COL_DISC), guide = "none") +
  scale_x_continuous(limits = c(-lim, lim)) + scale_y_continuous(limits = c(-lim, lim)) +
  labs(x = expression("Bulk RNA-seq " * log[2] * "FC"),
       y = expression("Hepatocyte scRNA " * log[2] * "FC")) +
  theme_masld_compact() + theme(aspect.ratio = 1)

message(sprintf(paste0("CAPTION (bulk<->scRNA hepatocyte concordance): per-gene bulk RNA-seq log2FC ",
  "(MASLD vs control) vs hepatocyte scRNA pseudobulk log2FC (MASLD vs Healthy), n=%d genes. ",
  "Dashed = perfect concordance (y=x). Gray = all genes (moderate concordance, Spearman rho=%.2f); ",
  "anchors (Fig 2/3/4 cast) coloured concordant [teal] vs discordant [magenta]. The scatter around ",
  "the diagonal is the honest per-gene view behind the single rho, and reflects composition-confound ",
  "+ cross-cohort differences (bulk and scRNA share no patients)."), nrow(dt), rho))

out <- file.path(PANEL_DIR, "bulk_sc_concordance_hep.pdf")
save_fig(p, out, width = fig_half_width, height = 3.0)
cat("[bulk<->sc hep] Saved:", out, "\n")
fwrite(dt[order(-abs(bulk_lfc))][, .(gene, bulk_lfc = round(bulk_lfc,3),
        sc_lfc = round(sc_lfc,3), anchor = !is.na(hl))],
       file.path(DATA_DIR, "bulk_sc_concordance_hep.csv"))
print(lab[, .(gene, bulk_lfc = round(bulk_lfc,2), sc_lfc = round(sc_lfc,2), hl)])
