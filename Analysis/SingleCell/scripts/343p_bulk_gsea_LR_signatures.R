#!/usr/bin/env Rscript
# 343p_bulk_gsea_LR_signatures.R
# ------------------------------------------------------------------------------
# Phase 7.5 — Bulk GSEA of stage-progressive LR-pair gene signatures.
#
# Score the stage-progressive LR ligand/receptor signatures (derived from sc-LMM)
# in the bulk Dream mega-analysis output as an orthogonal validation.
#
# Pre-registered threshold: NES > 1.5, p < 0.05 (per docs/archive/single_cell_pre_remediation_2026-05-22/scrna_pipeline_rigor_plan.md; superseded by docs/single_cell_analysis.md).
#
# Gene sets (per axis):
#   * `progressive_up`   = ligand∪receptor symbols of LR pairs with Estimate>0 AND padj_within_ct<0.10
#   * `progressive_down` = same with Estimate<0
#   * `dedup30_up`       = ligand∪receptor symbols from dedup-30 with Estimate>0
#   * `dedup30_down`     = same with Estimate<0
# Axes: SH (coarse), F-stage, pseudotime → 4 sets × 3 axes = 12 base sets;
# 12 sensitivity sets without HLA-A/HLA-B reported alongside.
#
# Ranking statistic: sign(logFC) × −log10(padj), per gene, on bulk Dream output.
# Bulk anchors: all-MASLD (dream_results_ashr.csv) + stage-specific (steatosis, sh, cirrhosis).
#
# Outputs:
#   Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/bulk_gsea_LR_signatures.tsv
#   figures/supplementary/stage_ccc/figS_bulk_gsea_LR.pdf
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(ggplot2)
  library(cowplot)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
                   "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))

STAGE_DIR <- file.path(BASE,
  "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory")
BULK_DIR  <- file.path(BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration")
FIG_DIR   <- file.path(BASE, "figures/supplementary/stage_ccc")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n", sep = " ")

set.seed(42)

# Pre-registered thresholds
NES_THRESH <- 1.5
P_THRESH   <- 0.05

# Significance cutoff for assembling LR gene sets (sc-LMM padj_within_ct)
LR_PADJ <- 0.10

# ──────────────────────────────────────────────────────────────────────────────
# 1. Load stage-LMM outputs and extract direction-signed LR gene sets
# ──────────────────────────────────────────────────────────────────────────────
log_msg("Loading stage-LMM LR results")

lmm_coarse     <- fread(file.path(STAGE_DIR, "stage_lr_lmm_coarse.tsv"))
lmm_fstage     <- fread(file.path(STAGE_DIR, "stage_lr_lmm_fstage.tsv"))
lmm_continuous <- fread(file.path(STAGE_DIR, "stage_lr_lmm_continuous.tsv"))
dedup30        <- fread(file.path(STAGE_DIR, "dedup_top30_independent_signals.tsv"))

log_msg(sprintf("coarse: %d rows | fstage: %d rows | pseudotime: %d rows | dedup30: %d rows",
                nrow(lmm_coarse), nrow(lmm_fstage), nrow(lmm_continuous), nrow(dedup30)))

# Filter terms per axis (only retain the canonical progression term)
coarse_sh <- lmm_coarse[term == "disease_stage_coarseSteatohepatitis"]
fstage    <- lmm_fstage[term == "F_stage_numeric"]
pseudo    <- lmm_continuous[term == "macrophage_pseudotime_mean"]

log_msg(sprintf("Per-axis filtered rows | coarse_SH=%d, F-stage=%d, pseudotime=%d",
                nrow(coarse_sh), nrow(fstage), nrow(pseudo)))

# Helper: build ligand∪receptor symbol set with direction sign + padj filter
build_lr_set <- function(dt, direction = c("up", "down"),
                          padj_col = "padj_within_ct", padj_thr = LR_PADJ) {
  direction <- match.arg(direction)
  sig <- dt[get(padj_col) < padj_thr]
  if (direction == "up") sig <- sig[Estimate > 0] else sig <- sig[Estimate < 0]
  syms <- unique(c(sig$ligand_complex, sig$receptor_complex))
  # Some receptors are complexes joined by "_" (rare) — split conservatively
  syms <- unique(unlist(strsplit(syms, "_", fixed = TRUE)))
  syms <- syms[!is.na(syms) & nchar(syms) > 0]
  syms
}

# dedup-30: padj already passed at sc-side, just use Estimate sign
build_dedup_set <- function(dt, direction = c("up", "down")) {
  direction <- match.arg(direction)
  sig <- if (direction == "up") dt[Estimate > 0] else dt[Estimate < 0]
  syms <- unique(c(sig$representative_ligand, sig$receptor))
  # Some receptors include subunits
  syms <- unique(unlist(strsplit(syms, "_", fixed = TRUE)))
  syms <- syms[!is.na(syms) & nchar(syms) > 0]
  syms
}

# dedup-30 axis tagging: we look up each lr_pair across the three axes to determine
# axis-specific direction. The dedup TSV does not have explicit axis, but the
# Estimate column is the axis chosen during ranking. For per-axis dedup gene sets,
# we re-anchor: filter lmm_<axis> for the dedup-30 lr_pairs and take direction there.
dedup_axis_set <- function(axis_dt, dedup_dt, direction) {
  pairs <- unique(dedup_dt$lr_pair)
  sub <- axis_dt[lr_pair %in% pairs & padj_within_ct < LR_PADJ]
  if (direction == "up") sub <- sub[Estimate > 0] else sub <- sub[Estimate < 0]
  syms <- unique(c(sub$ligand_complex, sub$receptor_complex))
  syms <- unique(unlist(strsplit(syms, "_", fixed = TRUE)))
  syms <- syms[!is.na(syms) & nchar(syms) > 0]
  syms
}

axis_tab <- list(
  SH        = coarse_sh,
  Fstage    = fstage,
  Pseudotime = pseudo
)

gene_sets <- list()
for (axis_name in names(axis_tab)) {
  dt <- axis_tab[[axis_name]]
  gene_sets[[sprintf("LR_%s_up", axis_name)]]      <- build_lr_set(dt, "up")
  gene_sets[[sprintf("LR_%s_down", axis_name)]]    <- build_lr_set(dt, "down")
  gene_sets[[sprintf("dedup30_%s_up", axis_name)]] <- dedup_axis_set(dt, dedup30, "up")
  gene_sets[[sprintf("dedup30_%s_down", axis_name)]] <- dedup_axis_set(dt, dedup30, "down")
}

# HLA-removal sensitivity sets
hla_strip <- c("HLA-A", "HLA-B")
gene_sets_no_hla <- lapply(gene_sets, function(g) setdiff(g, hla_strip))
names(gene_sets_no_hla) <- paste0(names(gene_sets), "_noHLA")

all_sets <- c(gene_sets, gene_sets_no_hla)

log_msg("Gene set sizes:")
for (nm in names(all_sets)) {
  log_msg(sprintf("  %-35s n=%d", nm, length(all_sets[[nm]])))
}

# ──────────────────────────────────────────────────────────────────────────────
# 2. Build ranking vectors from bulk anchors
# ──────────────────────────────────────────────────────────────────────────────
make_ranks <- function(dt) {
  # Use signed -log10(padj). Tiny padj clamp avoids Inf.
  dt[, padj := pmax(padj, 1e-300)]
  dt[, signed_log10p := sign(logFC) * -log10(padj)]
  dt <- dt[!is.na(symbol) & nchar(symbol) > 0 & is.finite(signed_log10p)]
  # Collapse duplicate symbols by max |stat|
  dt[, abs_stat := abs(signed_log10p)]
  dt <- dt[order(-abs_stat)]
  dt <- dt[!duplicated(symbol)]
  r <- setNames(dt$signed_log10p, dt$symbol)
  sort(r, decreasing = TRUE)
}

log_msg("Loading bulk anchors")
bulk_files <- list(
  AllMASLD   = file.path(BULK_DIR, "dream_results_ashr.csv"),
  Steatosis  = file.path(BULK_DIR, "dream_results_stage_steatosis.csv"),
  SH         = file.path(BULK_DIR, "dream_results_stage_sh.csv"),
  Cirrhosis  = file.path(BULK_DIR, "dream_results_stage_cirrhosis.csv")
)

ranks_by_anchor <- lapply(bulk_files, function(p) {
  d <- fread(p)
  # Stage-specific CSV uses padj column; all-MASLD ashr uses padj column too.
  if (!"padj" %in% names(d) && "P.Value" %in% names(d)) {
    # fall-back: use P.Value
    d[, padj := P.Value]
  }
  make_ranks(d)
})

for (nm in names(ranks_by_anchor)) {
  log_msg(sprintf("Anchor %-10s n=%d genes ranked", nm, length(ranks_by_anchor[[nm]])))
}

# ──────────────────────────────────────────────────────────────────────────────
# 3. Run fgsea per gene set × anchor
# ──────────────────────────────────────────────────────────────────────────────
log_msg("Running fgsea")

gsea_long <- list()
for (anchor_nm in names(ranks_by_anchor)) {
  r <- ranks_by_anchor[[anchor_nm]]
  fg <- fgseaMultilevel(
    pathways   = all_sets,
    stats      = r,
    minSize    = 3,
    maxSize    = 5000,
    nPermSimple = 1000,
    eps        = 0
  )
  fg <- as.data.table(fg)
  fg[, bulk_anchor := anchor_nm]
  gsea_long[[anchor_nm]] <- fg
}
gsea_dt <- rbindlist(gsea_long, fill = TRUE)

# Annotate axis / direction / set type / HLA flag
gsea_dt[, axis      := fifelse(grepl("_SH_", pathway), "SH",
                       fifelse(grepl("_Fstage_", pathway), "F_stage",
                       fifelse(grepl("_Pseudotime_", pathway), "Pseudotime", NA_character_)))]
gsea_dt[, direction := fifelse(grepl("_up", pathway), "up",
                       fifelse(grepl("_down", pathway), "down", NA_character_))]
gsea_dt[, set_type  := fifelse(grepl("^LR_", pathway), "LR",
                       fifelse(grepl("^dedup30_", pathway), "dedup30", NA_character_))]
gsea_dt[, hla_strip := grepl("_noHLA$", pathway)]

# BH correction *across the 12 fgsea tests* of primary gene sets (LR, full, no HLA strip)
# per the constraint. We compute padj_BH_primary separately for the 12 primary
# (LR × {up,down} × {SH,Fstage,Pseudotime}) tests in the AllMASLD anchor.
primary_mask <- with(gsea_dt,
                     set_type == "LR" & !hla_strip & bulk_anchor == "AllMASLD")
gsea_dt[, padj_BH_primary12 := NA_real_]
gsea_dt[primary_mask,
        padj_BH_primary12 := p.adjust(pval, method = "BH")]

# Verdict for matched-direction success
gsea_dt[, matched_success := (NES > NES_THRESH & pval < P_THRESH)]

# leadingEdge to scalar
gsea_dt[, leadingEdge_str := vapply(leadingEdge,
                                    function(x) paste(x, collapse = ";"),
                                    character(1))]

out_cols <- c("axis", "direction", "set_type", "hla_strip", "bulk_anchor",
              "pathway", "size", "ES", "NES", "pval", "padj",
              "padj_BH_primary12", "matched_success", "leadingEdge_str")
gsea_out <- gsea_dt[, ..out_cols]
setnames(gsea_out, c("padj", "size"), c("padj_fgsea", "gene_set_size"))
setorder(gsea_out, bulk_anchor, axis, set_type, direction, hla_strip)

fwrite(gsea_out,
       file.path(STAGE_DIR, "bulk_gsea_LR_signatures.tsv"),
       sep = "\t")
log_msg("Wrote bulk_gsea_LR_signatures.tsv")

# ──────────────────────────────────────────────────────────────────────────────
# 4. Per-axis verdict summary
# ──────────────────────────────────────────────────────────────────────────────
# Matched-direction expectation:
#   progressive_up   in disease-up bulk → positive NES (significant)
#   progressive_down in disease-down bulk → negative NES (significant)
# Bulk is run as disease-vs-control → high signed_log10p = disease-up gene.
# A positive NES on progressive_up = matched. A positive NES on progressive_down
# in a disease-up rank is NOT matched (we expect NEGATIVE NES on progressive_down).
matched_dt <- gsea_dt[set_type == "LR" & !hla_strip & bulk_anchor == "AllMASLD"]
matched_dt[, expected_sign := fifelse(direction == "up", 1, -1)]
matched_dt[, matched_axis := (sign(NES) == expected_sign &
                              abs(NES) > NES_THRESH &
                              pval < P_THRESH)]

log_msg("Per-axis matched verdict (AllMASLD anchor, primary LR sets only):")
for (ax in unique(matched_dt$axis)) {
  for (dir in c("up", "down")) {
    row <- matched_dt[axis == ax & direction == dir]
    if (nrow(row) == 0) next
    log_msg(sprintf(
      "  %-12s | dir=%-4s | NES=%6.3f | p=%.3g | padj_BH=%.3g | size=%d | matched=%s",
      ax, dir, row$NES, row$pval, row$padj_BH_primary12, row$size,
      row$matched_axis))
  }
}

# Stage-stratified matched verdict (cirrhosis is the strongest anchor for F-stage axis)
log_msg("Stage-stratified anchors (LR sets, primary, no HLA strip):")
stage_dt <- gsea_dt[set_type == "LR" & !hla_strip & bulk_anchor != "AllMASLD"]
for (ax in unique(stage_dt$axis)) {
  for (anchor in c("Steatosis", "SH", "Cirrhosis")) {
    for (dir in c("up", "down")) {
      row <- stage_dt[axis == ax & direction == dir & bulk_anchor == anchor]
      if (nrow(row) == 0) next
      exp_sign <- if (dir == "up") 1 else -1
      matched <- (sign(row$NES) == exp_sign &
                  abs(row$NES) > NES_THRESH &
                  row$pval < P_THRESH)
      log_msg(sprintf(
        "  axis=%-10s anchor=%-10s dir=%-4s NES=%6.3f p=%.3g size=%d matched=%s",
        ax, anchor, dir, row$NES, row$pval, row$size, matched))
    }
  }
}

# ──────────────────────────────────────────────────────────────────────────────
# 5. Multi-panel figure (compact, ≤8×8 in)
# ──────────────────────────────────────────────────────────────────────────────
log_msg("Building figS_bulk_gsea_LR.pdf")

# Panel A: Heatmap of NES across axis × bulk_anchor × direction (LR primary only)
heat_dt <- gsea_dt[set_type == "LR" & !hla_strip]
heat_dt[, axis_dir := sprintf("%s\n%s", axis, direction)]
heat_dt[, signif_lbl := fifelse(pval < 0.001, "***",
                       fifelse(pval < 0.01,  "**",
                       fifelse(pval < 0.05,  "*",  "")))]
heat_dt[, bulk_anchor_f := factor(bulk_anchor,
                                  levels = c("AllMASLD", "Steatosis", "SH", "Cirrhosis"))]
heat_dt[, axis_dir_f := factor(axis_dir,
        levels = c("SH\nup", "SH\ndown",
                   "F_stage\nup", "F_stage\ndown",
                   "Pseudotime\nup", "Pseudotime\ndown"))]

pA <- ggplot(heat_dt,
             aes(x = bulk_anchor_f, y = axis_dir_f, fill = NES)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f%s", NES, signif_lbl)),
            size = 2.8, color = "black") +
  scale_fill_gradient2(low = masld_colors$control,
                       mid = "white",
                       high = masld_colors$mash,
                       midpoint = 0,
                       limits = c(-max(abs(heat_dt$NES), na.rm = TRUE),
                                   max(abs(heat_dt$NES), na.rm = TRUE)),
                       na.value = "grey90",
                       name = "NES") +
  labs(x = "Bulk anchor (disease vs control)",
       y = "Axis × direction\n(scRNA LR signature)",
       title = "A. Bulk GSEA of stage-progressive LR signatures",
       subtitle = "* p<0.05; ** p<0.01; *** p<0.001 (uncorrected)") +
  theme_minimal(base_size = 8) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        plot.title = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 7, color = "grey40"),
        legend.position = "right")

# Panel B: NES barplot for AllMASLD anchor only, faceted by axis
bar_dt <- heat_dt[bulk_anchor == "AllMASLD"]
bar_dt[, axis_f := factor(axis, levels = c("SH", "F_stage", "Pseudotime"))]
bar_dt[, dir_f  := factor(direction, levels = c("up", "down"))]
pB <- ggplot(bar_dt,
             aes(x = dir_f, y = NES,
                 fill = NES > 0)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = c(-NES_THRESH, NES_THRESH),
             linetype = "dashed", color = "grey40", linewidth = 0.3) +
  geom_text(aes(label = sprintf("p=%.2g\nn=%d", pval, size)),
            vjust = -0.3, size = 2.4) +
  scale_fill_manual(values = c(`TRUE` = masld_colors$mash,
                               `FALSE` = masld_colors$control),
                    guide = "none") +
  facet_wrap(~ axis_f, nrow = 1) +
  labs(x = "Signature direction",
       y = "NES (AllMASLD)",
       title = "B. AllMASLD anchor — matched direction expectation") +
  theme_minimal(base_size = 8) +
  theme(strip.text = element_text(face = "bold"),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(face = "bold", size = 9),
        axis.text = element_text(color = "black"))

# Panel C: Enrichment plot for the strongest matched signal (NES highest, dir=up)
top_row <- heat_dt[direction == "up" & bulk_anchor == "AllMASLD"][order(-NES)][1]
gs_for_plot <- all_sets[[top_row$pathway]]
top_rank <- ranks_by_anchor[["AllMASLD"]]
pC <- plotEnrichment(gs_for_plot, top_rank) +
  labs(title = sprintf("C. Enrichment: %s (NES=%.2f, p=%.2g)",
                       top_row$pathway, top_row$NES, top_row$pval)) +
  theme_minimal(base_size = 8) +
  theme(plot.title = element_text(face = "bold", size = 9),
        axis.text = element_text(color = "black"))

# Panel D: Sensitivity — NES with/without HLA-A,B (primary anchor)
sens_dt <- gsea_dt[set_type == "LR" & bulk_anchor == "AllMASLD"]
sens_dt[, hla_lbl := fifelse(hla_strip, "no HLA-A,B", "with HLA-A,B")]
sens_dt[, axis_dir_f := factor(sprintf("%s_%s", axis, direction),
        levels = c("SH_up", "SH_down",
                   "F_stage_up", "F_stage_down",
                   "Pseudotime_up", "Pseudotime_down"))]
pD <- ggplot(sens_dt,
             aes(x = axis_dir_f, y = NES,
                 group = hla_lbl, color = hla_lbl)) +
  geom_line(linewidth = 0.4, color = "grey50") +
  geom_point(size = 2.2) +
  geom_hline(yintercept = c(-NES_THRESH, 0, NES_THRESH),
             linetype = c("dashed", "solid", "dashed"),
             color = "grey50", linewidth = 0.3) +
  scale_color_manual(values = c(`with HLA-A,B` = masld_colors$mash,
                                `no HLA-A,B` = masld_colors$control),
                     name = NULL) +
  labs(x = NULL, y = "NES",
       title = "D. HLA-A/-B strip sensitivity (AllMASLD)") +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        plot.title = element_text(face = "bold", size = 9),
        legend.position = "top",
        legend.text = element_text(size = 7))

# Assemble: 2×2 grid, 8×8 in
fig <- plot_grid(pA, pB, pC, pD,
                 ncol = 2, rel_heights = c(1.05, 0.95),
                 align = "hv", axis = "tblr")

fig_path <- file.path(FIG_DIR, "figS_bulk_gsea_LR.pdf")
ggsave(fig_path, fig, width = 8, height = 8, units = "in", device = cairo_pdf)
log_msg(sprintf("Wrote %s", fig_path))

# ──────────────────────────────────────────────────────────────────────────────
# 6. Verdict text summary printed at end of log
# ──────────────────────────────────────────────────────────────────────────────
cat("\n", strrep("=", 78), "\n", sep = "")
cat("VERDICTS (LR primary sets, no HLA strip)\n")
cat(strrep("=", 78), "\n", sep = "")
for (ax in c("SH", "F_stage", "Pseudotime")) {
  for (anchor in c("AllMASLD", "Steatosis", "SH", "Cirrhosis")) {
    rup   <- gsea_dt[set_type == "LR" & !hla_strip & axis == ax &
                     direction == "up" & bulk_anchor == anchor]
    rdown <- gsea_dt[set_type == "LR" & !hla_strip & axis == ax &
                     direction == "down" & bulk_anchor == anchor]
    if (nrow(rup) == 0 || nrow(rdown) == 0) next
    up_ok   <- rup$NES > NES_THRESH & rup$pval < P_THRESH
    down_ok <- rdown$NES < -NES_THRESH & rdown$pval < P_THRESH
    ok_lbl  <- if (up_ok && down_ok) "BOTH PASS"
               else if (up_ok)       "UP PASS"
               else if (down_ok)     "DOWN PASS"
               else                  "FAIL"
    cat(sprintf(
      "axis=%-10s anchor=%-10s NES_up=%6.3f (p=%.2g) NES_down=%6.3f (p=%.2g) → %s\n",
      ax, anchor, rup$NES, rup$pval, rdown$NES, rdown$pval, ok_lbl))
  }
}
cat(strrep("=", 78), "\n", sep = "")

log_msg("343p complete.")
