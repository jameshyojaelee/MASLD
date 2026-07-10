#!/usr/bin/env Rscript
# KEY MESSAGE: Four independent modalities converge on ONE regulatory logic across
# fibrosis stages, and colour groups them by PATHWAY so the same biology seen in two
# modalities shares one hue: the metabolic-identity program (HNF4A TF activity +
# fatty-acid/peroxisomal (FAO) Hotspot module + THRB/RORA expression) is lost [BLUE], the
# fibrogenic/disease-driver program (SMAD3 TF activity + ductular-injury (BICC1) Hotspot
# module) is gained [ORANGE],
# NF-kB inflammation (RELA) is gained [MAGENTA], and T2D/Wnt (TCF7L2) is gained [GOLD].
# Each line is nominated by an orthogonal modality (GWAS motif-disruption, COLOC,
# SCENIC+ regulon, scRNA Hotspot), named in the caption (not on-panel).
# ============================================================================
# tf_convergence_vsF0.R  — Fig 3I (cross-modal convergence, vs F0)
#
# Three line-tracks sharing the fig3e F0..F4 (Kleiner) stage axis. Colour = PATHWAY
# (shared across tracks); the nominating modality is named in the caption.
#   Track 1 (TF activity)      decoupleR norm_wmean z on F{x}-vs-F0 dream t-stats
#                              (DoRothEA A/B/C): HNF4A(down) / RELA,SMAD3,TCF7L2(up).
#   Track 2 (module signature) mean shrunk log2FC vs F0 of scRNA hepatocyte Hotspot
#                              module genes (bulk_replicated), scored on bulk:
#                              Fatty-acid/peroxisomal FAO (Hep-19, down) /
#                              Ductular injury BICC1 (Hep-20, up).
#                              scRNA F-stage is inference-limited, so the module GENE
#                              SIGNATURES are scored on bulk (real Kleiner F-stage),
#                              NOT the cross-cohort scRNA F-stage axis.
#   Track 3 (expression)       shrunk log2FC vs F0 for THRB + RORA (COLOC anchors with
#                              no DoRothEA regulon; THRB = resmetirom target).
#
# Output: FIG2_DIR/panels/fig3i_tf_convergence_vsF0.pdf   (FIG2_DIR = fig3_RNAseq;
#   Fig 3I in the A-I layout; shares the Kleiner F0-F4 x-axis with the
#   progression cascade (now 3E))
# Cache : RNA-seq/results/stratified_causal/stage_vs_f0_tf_activity.csv
# ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggrepel); library(patchwork)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

PANEL_DIR   <- file.path(FIG2_DIR, "panels")
DATA_DIR    <- file.path(PANEL_DIR, "data")
OUT_PDF     <- file.path(PANEL_DIR, "fig3i_tf_convergence_vsF0.pdf")
ACT_CSV     <- file.path(BASE, "RNA-seq/results/stratified_causal/stage_vs_f0_tf_activity.csv")
DREAM_STAGE <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/disease_signatures/fibrosis_stage_dream.csv")
SYMB_FILE   <- file.path(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv")
MODULE_F    <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2/hotspot_modules/hepatocytes/module_genes.tsv")

STAGES   <- c("F0","F1","F2","F3","F4")
PANEL_WIDTH_IN  <- 3.04
PANEL_HEIGHT_IN <- 1.99
# Modest right expansion reserves a compact gutter for the short end-of-line labels
# (bare names; evidence provenance moved to the caption). Smaller than before so F0..F4
# span most of the width and roughly register with the 3E cascade (3E = 0.04).
X_EXPAND <- ggplot2::expansion(mult = c(0.03, 0.27))

# ---- Pathway palette (shared across tracks) --------------------------------
COL_IDENTITY <- "#1565C0"   # metabolic identity  (HNF4A TF + HNF4A module)
COL_THRB     <- "#0D47A1"   # metabolic identity  (dark blue)
COL_RORA     <- "#1E88E5"   # metabolic identity  (medium blue)
COL_TGFB     <- "#E64A19"   # fibrogenic/disease-driver gain (SMAD3/TGF-beta TF + ductular-injury BICC1 module)
COL_NFKB     <- "#C9265E"   # NF-kB inflammation
COL_WNT      <- "#F6A800"   # T2D / Wnt

symbol_map <- function() {
  m <- unique(fread(SYMB_FILE, select = c("gene","symbol"))[!is.na(symbol) & symbol != ""])
  m[, gene_nov := sub("\\..*$", "", gene)]; m
}
add_xpos <- function(dt) { dt[, stage := factor(stage, levels = STAGES)]
  dt[, xpos := as.numeric(stage) - 1]; dt[] }
lab_repel <- function(d) geom_text_repel(data = d, aes(label = lab), hjust = 0,
  nudge_x = 0.08, size = GEOM_TEXT_6PT, fontface = "italic", color = "black",
  lineheight = 0.85, segment.size = 0.2, segment.color = "grey70",
  min.segment.length = 0, box.padding = 0.12, xlim = c(4.08, 5.05),
  force = 3, max.overlaps = Inf, seed = 42, show.legend = FALSE)
  # 2026-07-09: dropped direction="y" (was forcing all end-labels to a single
  # vertical stack with no x-freedom) -- at the 3.04x1.99in target height,
  # SMAD3/TCF7L2 (close F4 values) had no room to separate and collided.
  # Unrestricted x/y repel + a wider xlim + higher force lets close-value
  # labels stagger horizontally instead of stacking on top of each other.

# ---------------------------------------------------------------------------
# vs-F0 TF activity (cache-guarded; IDENTICAL method to figS_transition_tf_activity.R)
# ---------------------------------------------------------------------------
if (!file.exists(ACT_CSV)) {
  message("Cache miss -> computing vs-F0 TF activity (decoupleR + DoRothEA)...")
  suppressPackageStartupMessages({ library(decoupleR); library(dorothea) })
  d0 <- fread(DREAM_STAGE); mm <- symbol_map()
  d0[, gene_nov := sub("\\..*$", "", gene)]
  d0 <- merge(d0, mm[, .(gene_nov, symbol)], by = "gene_nov", all.x = TRUE)
  d0 <- d0[!is.na(symbol) & symbol != ""][order(-abs(t))][!duplicated(paste0(symbol, "_", contrast))]
  reg <- as.data.frame(dorothea_hs[dorothea_hs$confidence %in% c("A","B","C"), c("tf","target","mor")])
  res <- rbindlist(lapply(c("F1_vs_F0","F2_vs_F0","F3_vs_F0","F4_vs_F0"), function(cc) {
    sub <- d0[contrast == cc]
    r <- as.data.table(run_wmean(matrix(setNames(sub$t, sub$symbol), ncol = 1,
                                        dimnames = list(sub$symbol, cc)),
                                 net = reg, .source = "tf", .target = "target",
                                 .mor = "mor", times = 1000, minsize = 5))
    r <- r[statistic == "norm_wmean"]; r[, contrast := cc]; r
  }))
  setnames(res, "source", "tf")
  dir.create(dirname(ACT_CSV), recursive = TRUE, showWarnings = FALSE)
  fwrite(res[, .(tf, contrast, score, p_value, statistic)], ACT_CSV)
}

# Shared bulk stage DE (symbol-mapped) for tracks 2 & 3
d <- fread(DREAM_STAGE); m <- symbol_map()
d[, gene_nov := sub("\\..*$", "", gene)]
d <- merge(d, m[, .(gene_nov, symbol)], by = "gene_nov", all.x = TRUE)
d[, eff := ifelse(is.na(shrunk_logFC), logFC, shrunk_logFC)]

# ---------------------------------------------------------------------------
# TRACK 1 — TF activity vs F0 (decoupleR z)
# ---------------------------------------------------------------------------
CONV <- c("HNF4A","RELA","SMAD3","TCF7L2")
act <- fread(ACT_CSV)[tf %in% CONV]
act[, stage := sub("_vs_F0", "", contrast)]
A <- rbind(act[, .(tf, stage, score, p_value)],
           data.table(tf = CONV, stage = "F0", score = 0, p_value = NA_real_), fill = TRUE)
A[, tf := factor(tf, levels = CONV)]
A[, lab := as.character(tf)]   # bare TF name; nominating modality named in the caption
A <- add_xpos(A)
act_pal <- c(HNF4A = COL_IDENTITY, RELA = COL_NFKB, SMAD3 = COL_TGFB, TCF7L2 = COL_WNT)

p1 <- ggplot(A, aes(xpos, score, color = tf, group = tf)) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_line(linewidth = 0.9) +
  geom_point(data = A[!is.na(p_value)], size = 1.45) +
  lab_repel(A[stage == "F4"]) +
  scale_color_manual(values = act_pal, guide = "none") +
  scale_x_continuous(breaks = 0:4, labels = STAGES, expand = X_EXPAND) +
  scale_y_continuous(name = "TF activity\nvs F0 (z)", breaks = seq(-6, 12, 3)) +
  labs(x = NULL) + coord_cartesian(clip = "off") +
  theme_masld_compact() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank(), plot.margin = margin(2, 34,0, 3))

# ---------------------------------------------------------------------------
# TRACK 2 — scRNA Hotspot module signatures scored on bulk (mean log2FC vs F0)
# ---------------------------------------------------------------------------
mg <- fread(MODULE_F)
mod_def <- list(`19` = "Fatty-acid / peroxisomal (FAO)", `20` = "Ductular injury (BICC1)")
mod_col <- c(`Fatty-acid / peroxisomal (FAO)` = COL_IDENTITY, `Ductular injury (BICC1)` = COL_TGFB)
M <- rbindlist(lapply(names(mod_def), function(mi) {
  gg <- mg[module == as.integer(mi), gene]
  s <- d[symbol %in% gg, .(eff = mean(eff, na.rm = TRUE)),
         by = .(stage = sub("_vs_F0", "", contrast))]
  s[, mod := mod_def[[mi]]]; s
}))
M <- rbind(M, data.table(mod = unlist(mod_def), stage = "F0", eff = 0), fill = TRUE)
M[, mod := factor(mod, levels = unlist(mod_def))]
mod_lab <- c(`Fatty-acid / peroxisomal (FAO)` = "Fatty-acid oxidation", `Ductular injury (BICC1)` = "Ductular injury")
M[, lab := mod_lab[as.character(mod)]]   # short module name; marker gene + provenance in the caption
M <- add_xpos(M)

p2 <- ggplot(M, aes(xpos, eff, color = mod, group = mod)) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_line(linewidth = 0.9) +
  geom_point(data = M[stage != "F0"], size = 1.45) +
  lab_repel(M[stage == "F4"]) +
  scale_color_manual(values = mod_col, guide = "none") +
  scale_x_continuous(breaks = 0:4, labels = STAGES, expand = X_EXPAND) +
  scale_y_continuous(name = "Hotspot\nlog2FC", breaks = c(-0.5, 0, 0.5, 1.0)) +
  labs(x = NULL) + coord_cartesian(clip = "off") +
  theme_masld_compact() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank(), plot.margin = margin(0, 34,0, 3))

# ---------------------------------------------------------------------------
# TRACK 3 — genetic-anchor gene expression vs F0 (THRB, RORA)
# ---------------------------------------------------------------------------
GENES <- c("THRB","RORA")
ex <- d[symbol %in% GENES, .(gene = symbol, stage = sub("_vs_F0", "", contrast), eff, padj)]
ex[, sig := padj < 0.05]
E <- rbind(ex, data.table(gene = GENES, stage = "F0", eff = 0, padj = NA_real_, sig = NA), fill = TRUE)
E[, gene := factor(gene, levels = GENES)]
E[, lab := as.character(gene)]   # bare gene name; COLOC provenance in the caption
E <- add_xpos(E)
gene_pal <- c(THRB = COL_THRB, RORA = COL_RORA)

p3 <- ggplot(E, aes(xpos, eff, color = gene, group = gene)) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_line(linewidth = 0.9) +
  geom_point(data = E[!is.na(padj)], aes(shape = sig), size = 1.55, stroke = 0.4, fill = "white") +
  lab_repel(E[stage == "F4"]) +
  scale_color_manual(values = gene_pal, guide = "none") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21),
                     labels = c(`TRUE` = "padj<0.05", `FALSE` = "n.s."),
                     name = NULL, na.translate = FALSE) +
  scale_x_continuous(breaks = 0:4, labels = STAGES, expand = X_EXPAND) +
  scale_y_continuous(name = "Expression\nlog2FC", expand = expansion(mult = c(0.16, 0.12))) +
  labs(x = "Fibrosis stage (Kleiner F0–F4)") + coord_cartesian(clip = "off") +
  theme_masld_compact() +
  theme(axis.text.x = element_text(size = 6, face = "plain"),
        legend.position = c(0.17, 0.30), legend.key.size = unit(0.18, "cm"),
        legend.text = element_text(size = 6, face = "plain"), legend.background = element_blank(),
        plot.margin = margin(0, 34,3, 3))

# ---------------------------------------------------------------------------
# Assemble
# ---------------------------------------------------------------------------
panel <- p1 / p2 / p3 +
  plot_layout(heights = c(1.25, 0.9, 0.9)) +
  plot_annotation(
    title = "Multi-modal convergence across fibrosis stages",
    theme = theme(plot.title = element_text(size = 6, face = "plain", hjust = 0)))

ggsave(OUT_PDF, panel, width = PANEL_WIDTH_IN, height = PANEL_HEIGHT_IN, units = "in", device = cairo_pdf)

dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
fwrite(dcast(act, tf ~ stage, value.var = "score"),           file.path(DATA_DIR, "fig3i_conv_activity.csv"))
fwrite(dcast(M[stage != "F0"], mod ~ stage, value.var = "eff"), file.path(DATA_DIR, "fig3i_conv_module.csv"))
fwrite(dcast(ex, gene ~ stage, value.var = "eff"),            file.path(DATA_DIR, "fig3i_conv_expression.csv"))

message(sprintf(paste0(
  "[Fig 3I caption] Multi-modal convergence across fibrosis stages (F0 anchored at 0). ",
  "Colour = PATHWAY, shared across tracks: metabolic identity (blue), fibrogenic/disease-driver ",
  "gain (orange; SMAD3/TGF-beta TF activity + ductular-injury BICC1 module), NF-kB inflammation ",
  "(magenta), T2D/Wnt (gold). Each line's nominating orthogonal ",
  "modality is named below (provenance moved off-panel for legibility). ",
  "ALL THREE tracks are scored on the bulk RNA-seq fibrosis-stage DE (real per-sample ",
  "Kleiner F0-F4); the single-cell data contributes ONLY the Hotspot module gene sets ",
  "(defined by unsupervised gene-gene co-expression autocorrelation, which uses no stage ",
  "labels), never a single-cell stage axis.\n",
  "(1) TF activity (decoupleR norm_wmean z, DoRothEA A/B/C): HNF4A (top GWAS motif-disruption ",
  "load, 9 variants) falls to %.1f; RELA (NF-kB; COLOC 0.76), SMAD3 (TGF-beta; SCENIC+ regulon) ",
  "and TCF7L2 (T2D/Wnt; top motif load) rise. All p<0.05.\n",
  "(2) scRNA hepatocyte Hotspot module signatures scored on bulk (mean shrunk log2FC vs F0; ",
  "bulk_replicated): the fatty-acid/peroxisomal (FAO) module (ACOX2/EHHADH/PECR/LIPC) falls to %.2f -- ",
  "this FAO/peroxisomal loss IS the hepatocyte metabolic-identity program (HNF4A/PPARA-governed) -- ",
  "while the ductular-injury (BICC1) module, a disease-driver gain, rises to +%.2f.\n",
  "(3) Expression of the two strongest genetic anchors lacking a DoRothEA regulon: THRB (COLOC ",
  "0.9999, resmetirom target) falls to %.2f (padj 9.9e-13), RORA (COLOC 0.998) to %.2f. Filled ",
  "points padj<0.05. scRNA F-stage is inference-limited, so module signatures are scored on bulk ",
  "(real Kleiner F-stage), not the cross-cohort scRNA F-stage axis."),
  A[tf=="HNF4A" & stage=="F4", score],
  M[mod=="Fatty-acid / peroxisomal (FAO)" & stage=="F4", eff],
  M[mod=="Ductular injury (BICC1)" & stage=="F4", eff],
  E[gene=="THRB" & stage=="F4", eff], E[gene=="RORA" & stage=="F4", eff]))
cat(sprintf("[saved] %s\n", OUT_PDF))
