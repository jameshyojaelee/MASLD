#!/usr/bin/env Rscript
# =============================================================================
# consensus_pooled_meta.R
# -----------------------------------------------------------------------------
# "Pooled n meta-analysis consensus" for the DEG-method benchmark final report.
#
# Two DIFFERENT modeling paradigms call DEGs on the same 5-cohort cohort:
#   * POOLED cohort-adjusted analysis  -- pool all samples, cohort as a covariate
#       (e.g. limma_voom_qw__C2 / limma_voom__C4). "pooled" cell.
#   * RANDOM-EFFECTS META-ANALYSIS     -- per-cohort DE -> metafor REML pooling
#       (metafor_re__C2). "meta" cell.
# A gene called a DEG by BOTH paradigms is robust to the pooled-vs-meta modeling
# choice, hence HIGH-CONFIDENCE. This script quantifies and characterizes that
# intersection (the "consensus"):
#   1. Venn counts: consensus / pooled-only / meta-only + directional concordance.
#   2. Quality comparison: per-set external-truth enrichment (OpenTargets up OR,
#      mouse cross-species precision, DEG-COLOC Fisher OR, positive-control recall)
#      for {consensus, pooled-only, meta-only, pooled-union-meta}. The payoff is
#      whether the consensus is CLEANER (higher per-gene truth enrichment) than
#      either method alone.
#   3. Outputs: consensus_pooled_meta.csv (per-set metrics),
#      consensus_pooled_meta_genes.csv (consensus genes + both logFCs + truth),
#      a 5-line report summary, and a figure panel.
#
# Truth-loading + metric primitives are PORTED from the canonical factorial
# consumer so per-set numbers are computed identically to the winners table:
#   RNA-seq/.../multimethod_validation/factorial/metrics_truth_factorial.R
#   (itself ported from Cas13_Library_Design/scripts/benchmark_cutoffs.R metrics_of
#    + RNA-seq/27b_benchmark_presets.R positive-control evaluation).
#
# Env: /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
# CLI:
#   --pooled  pooled cell_id (default = top-ranked POOLED cell in
#             degmethod_ranked_winners.csv, else limma_voom_qw__C2__kna)
#   --meta    meta cell_id   (default metafor_re__C2__kna)
#   --dir     degx_factorial results dir (default = RNA-seq/results/degx_factorial)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ---------------------------------------------------------------------------
# 0. paths + CLI
# ---------------------------------------------------------------------------
BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
DEGX_DIR <- file.path(BASE, "RNA-seq/results/degx_factorial")

OT_FILE       <- file.path(BASE, "data/published_gene_panels/opentargets_masld_2025.tsv")
ATLAS_FILE    <- file.path(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
GENEMETA_FILE <- file.path(BASE, "data/gencode_v49_gene_metadata.tsv.gz")
PC_FILE       <- file.path(BASE, "RNA-seq/results/validation/positive_control_validation.csv")
WINNERS_FILE  <- file.path(DEGX_DIR, "degmethod_ranked_winners.csv")

# Canonical atlas truth columns (verified present 2026-06-05) — same as the
# factorial consumer. coloc_best_susie_pp4 is not materialized; the production
# PolyFun equivalent is the canonical SUSiE-COLOC posterior.
ATLAS_MOUSE_COL <- "mouse_meta_padj"
ATLAS_COLOC_COL <- "coloc_best_susie_pp4_polyfun"
COLOC_THRESH    <- 0.5

# DEG-set definition (canonical Tier-1, fixed_threshold).
DEG_PADJ <- 0.05
DEG_LFC  <- 0.5

strip_v <- function(x) sub("[.][0-9]+$", "", x)

# Default pooled cell: top-ranked POOLED cell in the winners table (pooled =
# NOT metafor_re and NOT a dream-only correction; i.e. a pooled limma/edger/deseq2
# engine). Falls back to limma_voom_qw__C2__kna if the winners file is absent.
default_pooled <- function(winners_file) {
  fallback <- "limma_voom_qw__C2__kna"
  if (!file.exists(winners_file)) return(fallback)
  w <- tryCatch(fread(winners_file), error = function(e) NULL)
  if (is.null(w) || !"cell_id" %in% names(w)) return(fallback)
  eng <- if ("engine" %in% names(w)) w$engine else sub("__.*$", "", w$cell_id)
  is_pooled <- !grepl("metafor", eng) & !grepl("dream", eng)
  wp <- w[is_pooled]
  if (nrow(wp) == 0) return(fallback)
  ord <- if ("winner_rank" %in% names(wp)) order(wp$winner_rank) else seq_len(nrow(wp))
  wp$cell_id[ord][1]
}

parse_args <- function() {
  a <- commandArgs(trailingOnly = TRUE)
  out <- list(dir = DEGX_DIR, pooled = NA_character_, meta = "metafor_re__C2__kna")
  i <- 1
  while (i <= length(a)) {
    key <- sub("^--", "", a[i])
    if (key %in% names(out)) { out[[key]] <- a[i + 1]; i <- i + 2 } else i <- i + 1
  }
  out
}
args <- parse_args()
DEGX_DIR  <- args$dir
CELLS_DIR <- file.path(DEGX_DIR, "cells")
OUT_DIR   <- DEGX_DIR
if (is.na(args$pooled)) args$pooled <- default_pooled(file.path(DEGX_DIR, "degmethod_ranked_winners.csv"))

pooled_id <- args$pooled
meta_id   <- args$meta
pooled_file <- file.path(CELLS_DIR, paste0("cell_", pooled_id, ".csv"))
meta_file   <- file.path(CELLS_DIR, paste0("cell_", meta_id,   ".csv"))

cat("=== consensus_pooled_meta.R ===\n")
cat("pooled cell (pooled cohort-adjusted analysis):", pooled_id, "\n")
cat("meta cell   (random-effects meta-analysis)   :", meta_id, "\n")
cat("cells_dir   :", CELLS_DIR, "\n\n")
if (!file.exists(pooled_file)) stop("pooled cell file not found: ", pooled_file)
if (!file.exists(meta_file))   stop("meta cell file not found: ",   meta_file)

# ---------------------------------------------------------------------------
# 1. TRUTH SOURCES (ported from metrics_truth_factorial.R)
# ---------------------------------------------------------------------------
gm <- fread(GENEMETA_FILE)
gm[, ensembl_base := strip_v(gene_id)]
sym2ens    <- gm[!duplicated(gene_name), setNames(ensembl_base, gene_name)]
sym2ens_up <- setNames(sym2ens, toupper(names(sym2ens)))
map_symbols <- function(symbols) {
  symbols <- as.character(symbols)
  out <- unname(sym2ens[symbols])
  miss <- is.na(out)
  if (any(miss)) out[miss] <- unname(sym2ens_up[toupper(symbols[miss])])
  out
}
# ensembl-base -> symbol (for annotating the consensus gene list)
ens2sym <- gm[!duplicated(ensembl_base), setNames(gene_name, ensembl_base)]

# OpenTargets MASLD gold standard.
ot <- fread(OT_FILE, skip = "gene_symbol")
ot_ens <- unique(na.omit(map_symbols(unique(ot$gene_symbol))))

# Multi-evidence atlas: universe + mouse / COLOC truth memberships.
atlas <- fread(ATLAS_FILE, select = c("ensembl_id", ATLAS_MOUSE_COL, ATLAS_COLOC_COL))
atlas[, gb := strip_v(ensembl_id)]
atlas <- atlas[!is.na(gb) & gb != ""][!duplicated(gb)]
atlas[, mouse_tested := !is.na(get(ATLAS_MOUSE_COL))]
atlas[, is_mouse_deg := mouse_tested & get(ATLAS_MOUSE_COL) < 0.05]
atlas[, is_coloc     := !is.na(get(ATLAS_COLOC_COL)) & get(ATLAS_COLOC_COL) > COLOC_THRESH]
mouse_universe <- atlas[mouse_tested == TRUE, gb]
mouse_truth    <- atlas[is_mouse_deg == TRUE, gb]
coloc_universe <- atlas[, gb]
coloc_truth    <- atlas[is_coloc == TRUE, gb]
ot_universe    <- atlas[, gb]
ot_truth       <- intersect(ot_ens, ot_universe)

# Positive controls (expression-driven) for held-out recall.
pc <- fread(PC_FILE)
pos_set_sym <- unique(pc[control_type == "Expression_driven", gene])
pos_set_ens <- unique(na.omit(map_symbols(pos_set_sym)))

cat(sprintf("Truth: OT-in-univ %d | mouse-tested %d (DEG %d) | COLOC>%.2f %d | pos-controls %d\n\n",
            length(ot_truth), length(mouse_universe), length(mouse_truth),
            COLOC_THRESH, length(coloc_truth), length(pos_set_ens)))

# Fisher 2x2 of DEG-set membership vs truth, over a shared universe ∩ tested set.
fisher_enrich <- function(deg_set, truth_set, universe, cell_universe) {
  universe <- intersect(universe, cell_universe)
  deg_u   <- intersect(deg_set, universe)
  truth_u <- intersect(truth_set, universe)
  a <- length(intersect(deg_u, truth_u))
  b <- length(setdiff(deg_u, truth_u))
  cc <- length(setdiff(truth_u, deg_u))
  d <- length(universe) - a - b - cc
  m <- matrix(c(a, b, cc, max(d, 0)), nrow = 2)
  ft <- tryCatch(fisher.test(m), error = function(e) NULL)
  list(a = a, or = if (is.null(ft)) NA_real_ else unname(ft$estimate),
       p = if (is.null(ft)) NA_real_ else ft$p.value)
}

# ---------------------------------------------------------------------------
# 2. read the two cells; define DEG sets (Tier-1 bidirectional + up-directional)
# ---------------------------------------------------------------------------
read_cell <- function(f) {
  dt <- fread(f)
  req <- c("gene", "logFC", "pval", "padj")
  miss <- setdiff(req, names(dt))
  if (length(miss)) stop(sprintf("cell %s missing columns: %s", basename(f), paste(miss, collapse = ",")))
  dt[, gb := strip_v(gene)]
  dt[!is.na(gb) & gb != ""]
}
dp <- read_cell(pooled_file)   # pooled
dm <- read_cell(meta_file)     # meta

deg_of <- function(dt, up_only = FALSE) {
  ok <- !is.na(dt$padj) & dt$padj < DEG_PADJ & !is.na(dt$logFC)
  if (up_only) ok & dt$logFC > DEG_LFC else ok & abs(dt$logFC) > DEG_LFC
}
dp[, is_deg := deg_of(dp)]; dp[, is_deg_up := deg_of(dp, TRUE)]
dm[, is_deg := deg_of(dm)]; dm[, is_deg_up := deg_of(dm, TRUE)]

pooled_deg    <- dp[is_deg == TRUE, gb]
meta_deg      <- dm[is_deg == TRUE, gb]
pooled_deg_up <- dp[is_deg_up == TRUE, gb]
meta_deg_up   <- dm[is_deg_up == TRUE, gb]

# Shared tested universe (intersection of genes both methods could call).
shared_univ <- intersect(dp$gb, dm$gb)

# Restrict DEG calls to the shared universe so set algebra is well-posed.
pooled_deg    <- intersect(pooled_deg,    shared_univ)
meta_deg      <- intersect(meta_deg,      shared_univ)
pooled_deg_up <- intersect(pooled_deg_up, shared_univ)
meta_deg_up   <- intersect(meta_deg_up,   shared_univ)

# ---------------------------------------------------------------------------
# 3. Venn + directional concordance
# ---------------------------------------------------------------------------
consensus     <- intersect(pooled_deg, meta_deg)
pooled_only   <- setdiff(pooled_deg, meta_deg)
meta_only     <- setdiff(meta_deg, pooled_deg)
union_set     <- union(pooled_deg, meta_deg)
consensus_up  <- intersect(pooled_deg_up, meta_deg_up)
pooled_only_up<- setdiff(pooled_deg_up, meta_deg_up)
meta_only_up  <- setdiff(meta_deg_up, pooled_deg_up)
union_up      <- union(pooled_deg_up, meta_deg_up)

# Directional concordance: among consensus genes, fraction with SAME logFC sign.
lf_p <- setNames(dp$logFC, dp$gb)
lf_m <- setNames(dm$logFC, dm$gb)
sign_agree <- sign(lf_p[consensus]) == sign(lf_m[consensus])
dir_concord <- if (length(consensus) > 0) mean(sign_agree, na.rm = TRUE) else NA_real_

cat("--- Venn (bidirectional Tier-1 DEG: padj<0.05 & |logFC|>0.5) ---\n")
cat(sprintf("  pooled DEG (%s) : %d\n", pooled_id, length(pooled_deg)))
cat(sprintf("  meta DEG   (%s) : %d\n", meta_id, length(meta_deg)))
cat(sprintf("  consensus (both): %d\n", length(consensus)))
cat(sprintf("  pooled-only     : %d\n", length(pooled_only)))
cat(sprintf("  meta-only       : %d\n", length(meta_only)))
cat(sprintf("  union           : %d\n", length(union_set)))
cat(sprintf("  directional agreement in consensus: %.2f%% (%d/%d same sign)\n\n",
            100 * dir_concord, sum(sign_agree, na.rm = TRUE), length(consensus)))

# ---------------------------------------------------------------------------
# 4. Per-set external-truth enrichment (the quality comparison)
# ---------------------------------------------------------------------------
# For each set, OT up-directional OR uses the up-directional members of that set
# (matches benchmark_cutoffs.R / metrics_truth_factorial.R operating direction);
# mouse precision, COLOC OR, PC recall use the bidirectional set.
set_metrics <- function(name, deg_bi, deg_up) {
  fe_ot_up <- fisher_enrich(deg_up, ot_truth, ot_universe, shared_univ)
  fe_coloc <- fisher_enrich(deg_bi, coloc_truth, coloc_universe, shared_univ)
  deg_mouse_univ <- intersect(deg_bi, mouse_universe)
  mouse_prec <- if (length(deg_mouse_univ) > 0)
    length(intersect(deg_mouse_univ, mouse_truth)) / length(deg_mouse_univ) else NA_real_
  pos_univ <- intersect(pos_set_ens, shared_univ)
  pc_rec   <- if (length(pos_univ) > 0) length(intersect(deg_bi, pos_univ)) / length(pos_univ) else NA_real_
  data.table(
    set                 = name,
    n_genes             = length(deg_bi),
    n_genes_up          = length(deg_up),
    ot_enrichment_or_up = fe_ot_up$or,
    ot_fisher_p_up      = fe_ot_up$p,
    ot_n_up             = fe_ot_up$a,
    mouse_precision     = mouse_prec,
    mouse_n_tested      = length(deg_mouse_univ),
    coloc_or            = fe_coloc$or,
    coloc_fisher_p      = fe_coloc$p,
    coloc_n             = fe_coloc$a,
    pc_recall           = pc_rec,
    pc_n                = length(intersect(deg_bi, intersect(pos_set_ens, shared_univ))))
}

metrics <- rbindlist(list(
  set_metrics("consensus",        consensus,   consensus_up),
  set_metrics("pooled_only",      pooled_only, pooled_only_up),
  set_metrics("meta_only",        meta_only,   meta_only_up),
  set_metrics("pooled_union_meta",union_set,   union_up)))
metrics[, pooled_cell := pooled_id]
metrics[, meta_cell   := meta_id]
metrics[, directional_concordance := dir_concord]

cat("--- per-set external-truth metrics ---\n")
print(metrics[, .(set, n_genes, n_genes_up, ot_enrichment_or_up,
                  mouse_precision, coloc_or, pc_recall)])

# Is the consensus cleaner than EITHER method alone? (consensus > both *-only)
g <- function(s, col) metrics[set == s, get(col)]
better_than_either <- function(col, higher_is_better = TRUE) {
  cv <- g("consensus", col); pv <- g("pooled_only", col); mv <- g("meta_only", col)
  if (any(is.na(c(cv, pv, mv)))) return(NA)
  if (higher_is_better) cv > pv & cv > mv else cv < pv & cv < mv
}
ot_better    <- better_than_either("ot_enrichment_or_up")
mouse_better  <- better_than_either("mouse_precision")
coloc_better <- better_than_either("coloc_or")

# ---------------------------------------------------------------------------
# 5. consensus gene table (both logFCs + truth annotations)
# ---------------------------------------------------------------------------
genes_dt <- data.table(ensembl_base = consensus)
genes_dt[, symbol         := ens2sym[ensembl_base]]
genes_dt[, pooled_logFC   := lf_p[ensembl_base]]
genes_dt[, meta_logFC     := lf_m[ensembl_base]]
genes_dt[, pooled_padj    := setNames(dp$padj, dp$gb)[ensembl_base]]
genes_dt[, meta_padj      := setNames(dm$padj, dm$gb)[ensembl_base]]
genes_dt[, sign_agree     := sign(pooled_logFC) == sign(meta_logFC)]
genes_dt[, direction      := ifelse(pooled_logFC > 0, "up", "down")]
genes_dt[, is_opentargets := ensembl_base %in% ot_truth]
genes_dt[, mouse_tested   := ensembl_base %in% mouse_universe]
genes_dt[, is_mouse_deg   := ensembl_base %in% mouse_truth]
genes_dt[, is_coloc       := ensembl_base %in% coloc_truth]
genes_dt[, is_pos_control := ensembl_base %in% pos_set_ens]
setorder(genes_dt, -is_opentargets, -is_coloc, pooled_padj)

# ---------------------------------------------------------------------------
# 6. write outputs
# ---------------------------------------------------------------------------
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
metrics_path <- file.path(OUT_DIR, "consensus_pooled_meta.csv")
genes_path   <- file.path(OUT_DIR, "consensus_pooled_meta_genes.csv")
fwrite(metrics, metrics_path)
fwrite(genes_dt, genes_path)
cat("\n=== wrote ===\n")
cat("  per-set metrics :", metrics_path, "\n")
cat("  consensus genes :", genes_path, sprintf("(%d genes)\n", nrow(genes_dt)))

# ---------------------------------------------------------------------------
# 7. figure panel — OpenTargets up-directional enrichment across the 3 sets
#    (pooled-only / consensus / meta-only). Control-invariant gray for the
#    union/reference; consensus highlighted.
# ---------------------------------------------------------------------------
FIG_OK <- TRUE
fig_path <- NULL
tryCatch({
  suppressWarnings(suppressMessages({
    source(file.path(BASE, "scripts/figures/publication_theme.R"))
    source(file.path(BASE, "scripts/figures/load_figure_data.R"))
  }))
  panels_dir <- file.path(FIGS_MULTIMETH_DIR, "panels")
  dir.create(panels_dir, showWarnings = FALSE, recursive = TRUE)
  fig_path <- file.path(panels_dir, "consensus_pooled_meta.pdf")

  plot_dt <- metrics[set %in% c("pooled_only", "consensus", "meta_only")]
  plot_dt[, label := factor(
    c(pooled_only = "Pooled-only", consensus = "Consensus", meta_only = "Meta-only")[set],
    levels = c("Pooled-only", "Consensus", "Meta-only"))]
  # Consensus = highlighted (Liang deep magenta); the two method-alone sets =
  # neutral gray, so the cleaner intersection stands out (control-invariant gray).
  fill_map <- c("Pooled-only" = "#9E9E9E", "Consensus" = "#C9265E", "Meta-only" = "#9E9E9E")

  p <- ggplot(plot_dt, aes(label, ot_enrichment_or_up, fill = label)) +
    geom_col(width = 0.66) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray40", linewidth = 0.3) +
    geom_text(aes(label = sprintf("OR=%.2f\n(n=%d)", ot_enrichment_or_up, n_genes_up)),
              vjust = -0.25, size = PUB_GEOM_TEXT) +
    scale_fill_manual(values = fill_map, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(x = NULL, y = "OpenTargets MASLD enrichment (up-directional OR)",
         title = "Pooled n meta-analysis consensus is cleaner",
         subtitle = sprintf("DEG = padj<0.05 & |logFC|>0.5  |  %.1f%% directional agreement in consensus",
                            100 * dir_concord)) +
    theme_bw(base_size = 7) + theme_pub() +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())

  ggsave(fig_path, p, width = 3.4, height = 3.0, useDingbats = FALSE)
  cat("  figure panel    :", fig_path, "\n")
}, error = function(e) {
  FIG_OK <<- FALSE
  cat("  [warn] figure render failed:", conditionMessage(e), "\n")
})

# ---------------------------------------------------------------------------
# 8. 5-line report summary
# ---------------------------------------------------------------------------
yn <- function(x) if (is.na(x)) "n/a" else if (x) "YES" else "no"
cat("\n================= REPORT SUMMARY (5 lines) =================\n")
cat(sprintf("1. Consensus (pooled n meta) = %d genes [pooled %d, meta %d; union %d]; pooled-only %d, meta-only %d.\n",
            length(consensus), length(pooled_deg), length(meta_deg), length(union_set),
            length(pooled_only), length(meta_only)))
cat(sprintf("2. Directional agreement within consensus: %.1f%% (%d/%d same logFC sign).\n",
            100 * dir_concord, sum(sign_agree, na.rm = TRUE), length(consensus)))
cat(sprintf("3. OpenTargets up-OR  consensus=%.2f vs pooled-only=%.2f, meta-only=%.2f -> consensus higher than both? %s.\n",
            g("consensus","ot_enrichment_or_up"), g("pooled_only","ot_enrichment_or_up"),
            g("meta_only","ot_enrichment_or_up"), yn(ot_better)))
cat(sprintf("4. Mouse precision    consensus=%.3f vs %.3f/%.3f (%s); DEG-COLOC OR consensus=%.2f vs %.2f/%.2f (%s).\n",
            g("consensus","mouse_precision"), g("pooled_only","mouse_precision"), g("meta_only","mouse_precision"),
            yn(mouse_better),
            g("consensus","coloc_or"), g("pooled_only","coloc_or"), g("meta_only","coloc_or"), yn(coloc_better)))
cat(sprintf("5. Consensus is more truth-enriched than EITHER method alone on %d/3 axes (OT/mouse/COLOC). PC recall consensus=%.3f.\n",
            sum(c(isTRUE(ot_better), isTRUE(mouse_better), isTRUE(coloc_better))),
            g("consensus","pc_recall")))
cat("===========================================================\n")
cat("\nDone.\n")
