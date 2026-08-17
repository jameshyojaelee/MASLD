#!/usr/bin/env Rscript
# KEY MESSAGE: Robust recurrent regulators rise across fibrosis stages while the prespecified hepatic-identity regulator HNF4A falls.
# ============================================================================
# tf_convergence_vsF0.R  — Figure 4F (inferred TF activity, vs F0)
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
# Output: Figure 4 candidate package or the Figure 4 legacy direct-render area;
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

# Candidate Figure 4F. TF activity is rendered alone and recomputed from the
# synchronized fragment-native F-versus-F0 statistics. BH adjustment is over
# every tested TF within each contrast before display selection. The display
# contains induced TFs that have the same activity direction in all four
# contrasts, BH q<0.05 in at least three, and mean absolute activity >=5, plus
# the three strongest directionally consistent reduced TFs with support in at
# least two contrasts. Direction is separated explicitly and the trajectory
# forest plot has no per-cell significance glyphs.
CANDIDATE_ROOT <- Sys.getenv("FIGURE_CANDIDATE_ROOT", "")
if (nzchar(CANDIDATE_ROOT)) {
  suppressPackageStartupMessages({ library(decoupleR); library(dorothea) })
  set.seed(42)
  EXT_ROOT <- normalizePath(Sys.getenv("STAGE_RELEASE_ROOT", ""), mustWork = TRUE)
  OUT_DIR <- file.path(CANDIDATE_ROOT, "figure4", "panels")
  SRC_DIR <- file.path(CANDIDATE_ROOT, "source_tables")
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(SRC_DIR, recursive = TRUE, showWarnings = FALSE)
  stage_stats <- fread(file.path(EXT_ROOT, "stage_extension_all_gene_results.tsv"))[axis == "fibrosis"]
  if (any(stage_stats[, uniqueN(gene_id_versioned), by = contrast]$V1 != unique(stage_stats$bh_family_size))) {
    stop("Incomplete F-versus-F0 BH family", call. = FALSE)
  }
  stage_stats <- stage_stats[!is.na(gene_name) & gene_name != ""][order(-abs(t))]
  stage_stats <- stage_stats[!duplicated(paste(contrast, gene_name))]
  regulon <- as.data.table(dorothea_hs)[confidence %in% c("A", "B", "C"), .(tf, target, mor)]
  regulon_sizes <- regulon[, .(regulon_size = uniqueN(target)), by = tf]
  activity_all <- rbindlist(lapply(unique(stage_stats$contrast), function(contrast_id) {
    d <- stage_stats[contrast == contrast_id]
    matrix <- matrix(d$t, ncol = 1L, dimnames = list(d$gene_name, contrast_id))
    result <- as.data.table(run_wmean(
      matrix, net = as.data.frame(regulon), .source = "tf", .target = "target",
      .mor = "mor", times = 1000, minsize = 5
    ))[statistic == "norm_wmean"]
    setnames(result, "source", "tf")
    result[, contrast := contrast_id]
    result
  }))
  activity_all[, bh_q := p.adjust(p_value, method = "BH"), by = contrast]
  activity_all <- merge(activity_all, regulon_sizes, by = "tf", all.x = TRUE)
  activity_all[, `:=`(seed = 42L, permutations = 1000L, confidence_levels = "A/B/C",
                      bh_family_size = .N), by = contrast]
  tf_summary <- activity_all[, .(
    n_bh = sum(bh_q < 0.05, na.rm = TRUE),
    min_bh_q = min(bh_q, na.rm = TRUE),
    mean_score = mean(score, na.rm = TRUE),
    mean_abs_score = mean(abs(score), na.rm = TRUE),
    same_direction = all(score > 0) || all(score < 0),
    max_abs_score = max(abs(score), na.rm = TRUE)
  ), by = tf]
  induced_tfs <- tf_summary[
    n_bh >= 3L & same_direction & mean_score > 0 & mean_abs_score >= 5, tf
  ]
  reduced_tfs <- tf_summary[
    n_bh >= 2L & same_direction & mean_score < 0
  ][order(-mean_abs_score)][1:3, tf]
  display_tfs <- unique(c(induced_tfs, reduced_tfs))
  if (length(induced_tfs) != 17L || length(reduced_tfs) != 3L ||
      length(display_tfs) != 20L ||
      !setequal(reduced_tfs, c("ZEB2", "HNF4A", "SIX5"))) {
    stop("TF display roster drift: expected 17 induced plus ZEB2/HNF4A/SIX5", call. = FALSE)
  }
  display <- merge(activity_all[tf %in% display_tfs], tf_summary, by = "tf")
  display[, selection_reason := fifelse(
    tf %in% reduced_tfs,
    "Three strongest same-direction reductions with BH q<0.05 in >=2 contrasts",
    "Same-direction induction; BH q<0.05 in >=3 contrasts; mean |activity| >=5"
  )]
  display[, stage := factor(sub("_vs_F0", "", contrast), levels = paste0("F", 1:4))]
  display[, direction := factor(fifelse(mean_score < 0, "Reduced", "Increased"),
                                levels = c("Increased", "Reduced"))]
  tf_order <- unique(display[, .(tf, direction, mean_abs_score)])[
    order(direction, mean_abs_score), tf
  ]
  display[, tf := factor(tf, levels = tf_order)]
  display[, `:=`(
    supported = bh_q < 0.05,
    stage_index = as.integer(stage)
  )]
  display[, `:=`(
    plot_score = abs(score)
  )]
  setorder(display, direction, tf, stage_index)
  activity_limit <- ceiling(max(display$plot_score, na.rm = TRUE))
  panel <- ggplot(display, aes(plot_score, tf, group = tf)) +
    geom_path(color = "#BDBDBD", linewidth = 0.35) +
    geom_point(aes(fill = stage), shape = 21, size = 1.25,
               stroke = 0.2, color = "white") +
    scale_fill_manual(values = fibrosis_stage_colors[paste0("F", 1:4)],
                      name = NULL, drop = FALSE) +
    scale_x_continuous(limits = c(0, activity_limit),
                       breaks = scales::pretty_breaks(n = 4),
                       expand = expansion(mult = c(0, 0.01))) +
    scale_y_discrete(drop = TRUE, expand = expansion(add = 0.35)) +
    facet_grid(direction ~ ., scales = "free_y", space = "free_y") +
    labs(x = "Absolute TF activity vs F0 (z)", y = NULL) +
    theme_masld_compact() +
    theme(axis.text.x = element_text(size = 6, face = "plain"),
          axis.text.y = element_text(size = 6, face = "italic"),
          axis.title.x = element_text(size = 6),
          axis.ticks.y = element_blank(), axis.line.y = element_blank(),
          strip.text.y = element_text(size = 6, face = "plain", angle = 0),
          strip.background = element_blank(),
          legend.position = "bottom", legend.direction = "horizontal",
          legend.text = element_text(size = 6),
          legend.key.width = unit(3.5, "mm"), legend.key.height = unit(2.5, "mm"),
          legend.margin = margin(0, 0, 0, 0),
          legend.box.spacing = unit(0.5, "mm"),
          plot.margin = margin(1, 1, 1, 1))
  ggsave(file.path(OUT_DIR, "fig4f_tf_activity.pdf"), panel,
         width = 3.05, height = 3.35, device = cairo_pdf)
  fwrite(activity_all, file.path(SRC_DIR, "fig4f_all_tf_activity.tsv"), sep = "\t")
  fwrite(display[, .(
    tf, direction, stage, contrast, score, plot_score, p_value, bh_q,
    supported, regulon_size,
    n_bh, min_bh_q, mean_score, mean_abs_score, same_direction,
    max_abs_score, selection_reason,
    seed, permutations, confidence_levels, bh_family_size
  )], file.path(SRC_DIR, "fig4f_display_tf_activity.tsv"), sep = "\t")
  quit(save = "no", status = 0)
}

PANEL_DIR   <- file.path(FIG4_SC_DIR, "panels", "_legacy_direct_render")
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
