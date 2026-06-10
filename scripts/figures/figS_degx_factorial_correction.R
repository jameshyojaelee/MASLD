#!/usr/bin/env Rscript
# figS_degx_factorial_correction.R
# -----------------------------------------------------------------------------
# Methods-robustness supplementary figure (panels J1-J5) for the FACTORIAL
# DEG-method benchmark: (engine x correction x k_sv) cells scored under the
# FROZEN calibration-gated composite rule (preregistration.yaml /
# build_ranked_table.R). SIBLING of figS_degx_multimethod.R -- shares its design
# language (engine-FAMILY categorical colour axis, [0.04,0.06] calibration band,
# N x N / lollipop conventions, data-driven `rd()` skip) but answers a different
# question: not "which engine", but "of all (engine x correction x k_sv) cells,
# which survived the type-I gate and ranked best, and is that robust to the
# batch-correction and surrogate-variable knobs?".
#
# DESIGN LANGUAGE (matches the multimethod battery):
#   * Engine family (deseq2 / edger / limma / dream / metafor) = the categorical
#     colour axis -- NEVER one colour per raw cell_id.
#   * Calibration band [CAL_LO=0.04, CAL_HI=0.06] shaded grey as the PASS zone.
#   * Control / baseline / neutral = #9E9E9E (baseline correction C0, gated-out
#     cells, the "be" anchor line).
#   * Sorted horizontal lollipops, no per-point value labels (colour + axis carry
#     it). PDF only (cairo_pdf, useDingbats=FALSE). ASCII. CPU job.
#   * Data-driven: each panel is guarded by an input-exists check; render what is
#     present, message what is skipped. The real metric tables do not all exist
#     yet, so the script MUST run gracefully on partial inputs.
#
# INPUTS (frozen schema, RNA-seq/results/degx_factorial/ -- override with --dir):
#   degmethod_ranked_winners.csv   one row per cell (build_ranked_table.R)
#   degmethod_scoring_long.csv     cell_id x metric x value x family x role
#   metrics_loco.csv               held-out LOCO reproducibility (RANK_R1)
#   metrics_external.csv           OT / mouse / coloc (RANK_R2)
#   metrics_calibration_..__*.csv  GLOB + rbind (per-cell type-I)
#   manifest_disease_vs_control*.csv  GLOB (n_sig_lfc05, k_sv, correction_id)
#   cells/cell_<id>.csv            DEG sets (gene,logFC,SE,pval,padj) for J5
#
# Run: micromamba run -n rnaseq Rscript scripts/figures/figS_degx_factorial_correction.R
# DO NOT COMMIT outputs (repo no-auto-commit rule).
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
})

BASE <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
source(file.path(BASE, "scripts/figures/publication_theme.R"))
source(file.path(BASE, "scripts/figures/load_figure_data.R"))

# --- CLI: --dir <degx_factorial dir> ---------------------------------------
parse_args <- function() {
  a <- commandArgs(trailingOnly = TRUE); out <- list(dir = NA_character_)
  i <- 1
  while (i <= length(a)) {
    key <- sub("^--", "", a[i])
    if (key %in% names(out)) { out[[key]] <- a[i + 1]; i <- i + 2 } else i <- i + 1
  }
  out
}
args <- parse_args()
DEGX_DIR <- if (!is.na(args$dir)) args$dir else file.path(BASE, "RNA-seq/results/degx_factorial")
CELLS_DIR <- file.path(DEGX_DIR, "cells")
CONTRAST  <- "disease_vs_control"

# --- constants (FROZEN; match figS_degx_multimethod.R + build_ranked_table.R)
NOMINAL <- 0.05; CAL_LO <- 0.04; CAL_HI <- 0.06
DEG_PADJ <- 0.05; DEG_LFC <- 0.5                  # DEG = padj<0.05 & |logFC|>0.5
SV_CORRECTIONS <- c("C4","C5","C6g","C6s","C6r","C6hk","C9")  # latent-factor / SV-bearing

OUT <- file.path(FIGS_MULTIMETH_DIR, "panels")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

cat("=== figS_degx_factorial_correction.R (panels J1-J5) ===\n")
cat("dir :", DEGX_DIR, "\n")
cat("out :", OUT, "\n\n")

# --- engine FAMILY map + palette (categorical colour axis) -----------------
fam_engine <- function(e){ e <- as.character(e)
  ifelse(grepl("^deseq2", e), "DESeq2",
  ifelse(grepl("^edger",  e), "edgeR",
  ifelse(grepl("^limma",  e), "limma",
  ifelse(grepl("^dream",  e), "dream",
  ifelse(grepl("^metafor",e), "metafor", "other")))))}
fam_cols <- c(DESeq2="#4aa2c2", edgeR="#518dc9", limma="#9b75d6",
              dream="#40b499", metafor="#e37faf", other="#BDBDBD")
FAM_LEVELS <- c("DESeq2","edgeR","limma","dream","metafor","other")

# pretty correction labels (short)
corr_labels <- c(
  C0="C0 none", C1="C1 tech-cov", C2="C2 bio-cov", C3="C3 tech+bio",
  C4="C4 cohort-FE", C5="C5 cohort-RE", C6g="C6g SVA", C6s="C6s SVA-2step",
  C6r="C6r RUVr", C6hk="C6hk RUVg-HK", C7="C7 ComBat-seq",
  C8="C8 quantile", C9="C9 PEER-PC", C10="C10 paired")
pc <- function(x){ x <- as.character(x)
  ifelse(x %in% names(corr_labels), corr_labels[x], x) }

# pretty engine labels
eng_labels <- c(
  deseq2_wald="DESeq2 (Wald)", deseq2_lrt="DESeq2 (LRT)",
  edger_qlf="edgeR (QLF)", edger_qlf_robust="edgeR (QLF-rob)", edger_lrt="edgeR (LRT)",
  limma_voom="limma-voom", limma_voom_qw="limma-voom (QW)", limma_trend="limma-trend",
  dream="dream", metafor_re="metafor (RE)")
pe <- function(x){ x <- as.character(x)
  ifelse(x %in% names(eng_labels), eng_labels[x], x) }

# k_sv pretty
pk <- function(x){ x <- as.character(x); sub("^k", "", x) }   # kbe->be, k10->10, kna->na

# --- shared helpers --------------------------------------------------------
rd <- function(...) { f <- file.path(...); if (file.exists(f)) utils::read.csv(f, check.names=FALSE, stringsAsFactors=FALSE) else NULL }
green_sc <- colorRampPalette(c("#dae7c7","#7fb069","#193c1e"))(12)
written <- character(0); mark <- function(p){written<<-c(written,p); message("  wrote ", basename(p))}
skip <- function(panel, why) message(sprintf("[skip] %s -- %s", panel, why))

# Derive engine/correction_id/k_sv from cell_id when absent (format eng__corr__k)
split_cell <- function(cell_id){
  parts <- strsplit(as.character(cell_id), "__", fixed=TRUE)
  data.frame(
    engine        = vapply(parts, function(p) if(length(p)>=1) p[1] else NA_character_, ""),
    correction_id = vapply(parts, function(p) if(length(p)>=2) p[2] else NA_character_, ""),
    k_sv          = vapply(parts, function(p) if(length(p)>=3) sub("^k","",p[3]) else NA_character_, ""),
    stringsAsFactors = FALSE)
}
# Backfill engine/correction_id/k_sv if a frame is missing them.
backfill_ids <- function(df){
  if (!"cell_id" %in% names(df)) return(df)
  s <- split_cell(df$cell_id)
  if (!"engine"        %in% names(df) || all(is.na(df$engine)))        df$engine        <- s$engine
  if (!"correction_id" %in% names(df) || all(is.na(df$correction_id))) df$correction_id <- s$correction_id
  if (!"k_sv"          %in% names(df) || all(is.na(df$k_sv)))          df$k_sv          <- s$k_sv
  # normalize k_sv (strip leading k if present, e.g. "kbe" -> "be")
  df$k_sv <- sub("^k", "", as.character(df$k_sv))
  df
}

# GLOB + rbind the per-cell calibration files.
load_calibration <- function(){
  pat <- sprintf("^metrics_calibration_%s__.*\\.csv$", CONTRAST)
  cfs <- list.files(DEGX_DIR, pattern = pat, full.names = TRUE)
  if (length(cfs) == 0) return(NULL)
  parts <- lapply(cfs, function(f) tryCatch(utils::read.csv(f, check.names=FALSE, stringsAsFactors=FALSE),
                                            error = function(e) NULL))
  parts <- Filter(function(d) !is.null(d) && nrow(d) > 0, parts)
  if (length(parts) == 0) return(NULL)
  cal <- dplyr::bind_rows(parts)
  if ("cell_id" %in% names(cal)) cal <- cal[!duplicated(cal$cell_id), , drop=FALSE]
  cal
}
# GLOB + rbind the manifest shards.
load_manifest <- function(){
  mfs <- list.files(DEGX_DIR, pattern = "^manifest_disease_vs_control.*\\.csv$", full.names = TRUE)
  if (length(mfs) == 0) return(NULL)
  parts <- lapply(mfs, function(f) tryCatch(utils::read.csv(f, check.names=FALSE, stringsAsFactors=FALSE),
                                            error = function(e) NULL))
  parts <- Filter(function(d) !is.null(d) && nrow(d) > 0, parts)
  if (length(parts) == 0) return(NULL)
  # rbindlist(fill=TRUE) tolerates heterogeneously-typed manifest shards (coerces
  # to a common type); dplyr::bind_rows() errors on type mismatch. Matches the scorer.
  mf <- as.data.frame(data.table::rbindlist(parts, use.names = TRUE, fill = TRUE),
                      stringsAsFactors = FALSE)
  if ("cell_id" %in% names(mf)) {
    if ("status" %in% names(mf)) mf <- mf[order(mf$cell_id, mf$status != "ok"), , drop=FALSE]
    mf <- mf[!duplicated(mf$cell_id), , drop=FALSE]
  }
  mf
}

# ---- load the primary winners table (drives J1, J2; informs J4) -----------
W <- rd(DEGX_DIR, "degmethod_ranked_winners.csv")
if (!is.null(W)) {
  W <- backfill_ids(W)
  W$fam <- factor(fam_engine(W$engine), levels = FAM_LEVELS)
  # coerce logical-ish columns robustly
  asbool <- function(x) { if (is.logical(x)) return(x); tolower(as.character(x)) %in% c("true","t","1","yes") }
  for (cc in c("eligible","gate_pass","is_winner"))
    if (cc %in% names(W)) W[[cc]] <- asbool(W[[cc]])
}

# ===========================================================================
# J1. CALIBRATION GATE -- who survived, and of those, who ranked best.
#     x = perm_typeI_mean (PASS band [0.04,0.06] shaded);
#     y = composite_ranksum (LOWER = better -> y-axis reversed);
#     gated-out / ineligible cells greyed; survivors coloured by engine family;
#     winner annotated.
# ===========================================================================
if (is.null(W) || !"perm_typeI_mean" %in% names(W) || all(is.na(W$perm_typeI_mean))) {
  skip("J1 calibration gate", "degmethod_ranked_winners.csv absent or no perm_typeI_mean")
} else {
  message("[render] J1 calibration gate")
  d <- W
  d$gp <- d$gate_pass %in% TRUE
  d$el <- if ("eligible" %in% names(d)) d$eligible %in% TRUE else TRUE
  d$survivor <- d$gp & d$el
  # survivors get a composite_ranksum; cells without one (ungated/ineligible) are
  # plotted along a baseline row so the gate decision still reads.
  has_comp <- "composite_ranksum" %in% names(d) && any(is.finite(d$composite_ranksum))
  if (!has_comp) d$composite_ranksum <- NA_real_
  base_y <- if (any(is.finite(d$composite_ranksum))) max(d$composite_ranksum, na.rm=TRUE) * 1.08 else 1
  d$y <- ifelse(is.finite(d$composite_ranksum), d$composite_ranksum, base_y)

  surv <- d[d$survivor & is.finite(d$perm_typeI_mean), , drop=FALSE]
  rest <- d[!(d$survivor) & is.finite(d$perm_typeI_mean), , drop=FALSE]
  win  <- d[("is_winner" %in% names(d)) & (d$is_winner %in% TRUE), , drop=FALSE]

  xr <- range(c(d$perm_typeI_mean, CAL_LO, CAL_HI), na.rm=TRUE)
  pJ1 <- ggplot() +
    annotate("rect", xmin=CAL_LO, xmax=CAL_HI, ymin=-Inf, ymax=Inf, fill="#9E9E9E", alpha=0.14) +
    geom_vline(xintercept=NOMINAL, linetype="dotted", linewidth=0.25, color="grey45")
  if (nrow(rest))
    pJ1 <- pJ1 + geom_point(data=rest, aes(perm_typeI_mean, y),
                            color="#9E9E9E", size=1.3, alpha=0.6, shape=16)
  if (nrow(surv))
    pJ1 <- pJ1 + geom_point(data=surv, aes(perm_typeI_mean, y, color=fam), size=2.0)
  if (nrow(win)) {
    pJ1 <- pJ1 + geom_point(data=win, aes(perm_typeI_mean, y), shape=21,
                            color="black", fill=NA, size=3.4, stroke=0.6)
    wlab <- sprintf("winner: %s", win$cell_id[1])
    lab_layer <- if (requireNamespace("ggrepel", quietly=TRUE))
        ggrepel::geom_text_repel(data=win, aes(perm_typeI_mean, y, label=wlab),
          size=PUB_GEOM_TEXT, color="grey15", seed=42, min.segment.length=0,
          segment.size=0.2, box.padding=0.6)
      else geom_text(data=win, aes(perm_typeI_mean, y, label=wlab),
          size=PUB_GEOM_TEXT, color="grey15", vjust=-1)
    pJ1 <- pJ1 + lab_layer
  }
  pJ1 <- pJ1 +
    scale_color_manual(values=fam_cols, name="engine family", drop=TRUE) +
    scale_y_reverse() +
    coord_cartesian(xlim=xr) +
    labs(title="Calibration gate: who survives, and which survivor ranks best?",
         subtitle="Grey band [0.04,0.06] = type-I PASS zone; grey points = gated-out / ineligible cells; coloured = survivors; ring = winner",
         x="permutation type-I error (per cell)",
         y="composite rank-sum  (lower = better; top = best)") +
    theme_masld()+theme_pub()+theme(legend.position="right")
  save_fig(pJ1, file.path(OUT, "panelJ1_calibration_gate.pdf"),
           width=fig_col_width, height=fig_col_width*0.82)
  mark(file.path(OUT, "panelJ1_calibration_gate.pdf"))
}

# ===========================================================================
# J2. COMPOSITE RANK -- horizontal lollipops of gate-pass AND eligible cells,
#     sorted by composite_ranksum, coloured by engine family. R1_rank + R2_rank
#     shown as small adjacent ticks so the reader sees WHY a cell won.
# ===========================================================================
if (is.null(W) || !"composite_ranksum" %in% names(W) || all(!is.finite(W$composite_ranksum))) {
  skip("J2 composite rank", "degmethod_ranked_winners.csv absent or no composite_ranksum (rank pending)")
} else {
  message("[render] J2 composite rank")
  el <- if ("eligible" %in% names(W)) W$eligible %in% TRUE else TRUE
  gp <- W$gate_pass %in% TRUE
  d <- W[gp & el & is.finite(W$composite_ranksum), , drop=FALSE]
  if (nrow(d) == 0) {
    skip("J2 composite rank", "no gate-pass & eligible cells with a composite rank")
  } else {
    d <- d[order(-d$composite_ranksum), , drop=FALSE]   # best at top after factor
    d$lab <- factor(d$cell_id, levels = d$cell_id)
    have_r <- all(c("R1_rank","R2_rank") %in% names(d))
    pJ2 <- ggplot(d, aes(composite_ranksum, lab)) +
      geom_segment(aes(x=0, xend=composite_ranksum, y=lab, yend=lab),
                   color="grey80", linewidth=0.3) +
      geom_point(aes(color=fam), size=2.1)
    if (have_r) {
      r_long <- d %>%
        select(lab, R1_rank, R2_rank) %>%
        tidyr::pivot_longer(c(R1_rank, R2_rank), names_to="which", values_to="rk")
      pJ2 <- pJ2 +
        geom_point(data=r_long, aes(rk, lab, shape=which),
                   color="grey35", size=1.0, inherit.aes=FALSE) +
        scale_shape_manual(values=c(R1_rank=124, R2_rank=43),
                           labels=c(R1_rank="R1 (LOCO)", R2_rank="R2 (external)"),
                           name="component rank")
    }
    if (any(d$is_winner %in% TRUE)) {
      wlab <- d[d$is_winner %in% TRUE, , drop=FALSE]
      pJ2 <- pJ2 + geom_point(data=wlab, aes(composite_ranksum, lab),
                              shape=21, color="black", fill=NA, size=3.4, stroke=0.6)
    }
    pJ2 <- pJ2 +
      scale_color_manual(values=fam_cols, name="engine family", drop=TRUE) +
      scale_y_discrete(labels=function(x) x) +
      scale_x_continuous(expand=expansion(mult=c(0,0.06))) +
      labs(title="Composite ranking of gate-passing, eligible cells",
           subtitle="Dot = composite rank-sum (lower = better); grey ticks = R1 (held-out LOCO) and R2 (external) component ranks; ring = winner",
           x="composite rank-sum (lower = better)", y=NULL) +
      theme_masld()+theme_pub()+theme(legend.position="right")
    save_fig(pJ2, file.path(OUT, "panelJ2_composite_rank.pdf"),
             width=fig_col_width, height=max(3.0, 0.16*nrow(d)+1.2))
    mark(file.path(OUT, "panelJ2_composite_rank.pdf"))
  }
}

# ===========================================================================
# J3. CORRECTION x DISEASE AXIS -- held-out replication (loco_repro_scalar)
#     distribution PER correction_id, across engines. Shows whether the
#     batch-correction choice moves held-out performance (expected second-order).
#     Baseline correction C0 in #9E9E9E.
# ===========================================================================
{
  # Prefer the winners table (has loco_repro_scalar + ids); fall back to
  # metrics_loco.csv (loco_repro_scalar or loco_auc) joined to ids by cell_id.
  d3 <- NULL
  if (!is.null(W) && "loco_repro_scalar" %in% names(W) && any(is.finite(W$loco_repro_scalar))) {
    d3 <- W[, intersect(c("cell_id","engine","correction_id","k_sv","loco_repro_scalar"), names(W)), drop=FALSE]
    d3$metric_val <- d3$loco_repro_scalar; d3$metric_name <- "loco_repro_scalar"
  } else {
    loco <- rd(DEGX_DIR, "metrics_loco.csv")
    if (!is.null(loco) && "cell_id" %in% names(loco)) {
      vcol <- if ("loco_repro_scalar" %in% names(loco)) "loco_repro_scalar"
              else if ("loco_auc" %in% names(loco)) "loco_auc" else NA_character_
      if (!is.na(vcol)) {
        loco <- backfill_ids(loco)
        d3 <- loco[, intersect(c("cell_id","engine","correction_id","k_sv",vcol), names(loco)), drop=FALSE]
        d3$metric_val <- d3[[vcol]]; d3$metric_name <- vcol
      }
    }
  }
  if (is.null(d3) || all(!is.finite(d3$metric_val))) {
    skip("J3 correction x disease axis", "no loco_repro_scalar / loco_auc available")
  } else {
    message("[render] J3 correction x disease axis (", unique(d3$metric_name)[1], ")")
    d3 <- d3[is.finite(d3$metric_val), , drop=FALSE]
    d3$fam  <- factor(fam_engine(d3$engine), levels=FAM_LEVELS)
    # order corrections by median held-out performance; C0 always rendered grey
    ord <- d3 %>% group_by(correction_id) %>%
      summarise(m=median(metric_val, na.rm=TRUE), .groups="drop") %>%
      arrange(m) %>% pull(correction_id)
    d3$corr <- factor(pc(d3$correction_id), levels=pc(ord))
    is_c0 <- function(lv) grepl("^C0\\b", lv)
    box_fill <- ifelse(is_c0(levels(d3$corr)), "#9E9E9E", "white")
    pJ3 <- ggplot(d3, aes(corr, metric_val)) +
      geom_boxplot(aes(fill=corr), outlier.shape=NA, width=0.6, linewidth=0.3, color="grey45") +
      geom_jitter(aes(color=fam), width=0.14, height=0, size=1.4, alpha=0.9) +
      scale_fill_manual(values=setNames(box_fill, levels(d3$corr)), guide="none") +
      scale_color_manual(values=fam_cols, name="engine family", drop=TRUE) +
      labs(title="Does the batch-correction choice move held-out replication?",
           subtitle="Held-out LOCO reproducibility per correction (each point = one engine x k_sv cell); baseline C0 box in grey. Flat = correction is second-order.",
           x="correction", y=unique(d3$metric_name)[1]) +
      theme_masld()+theme_pub()+
      theme(axis.text.x=element_text(angle=40, hjust=1), legend.position="right")
    save_fig(pJ3, file.path(OUT, "panelJ3_correction_disease_axis.pdf"),
             width=fig_full_width, height=3.6)
    mark(file.path(OUT, "panelJ3_correction_disease_axis.pdf"))
  }
}

# ===========================================================================
# J4. SV SWEEP -- for SVA/RUV/latent-factor corrections, scan k_sv in
#     {be,2,5,10,20}: y = three small-multiple lines {perm_typeI,
#     loco_repro_scalar, ot_enrichment_or_up}, faceted/coloured by correction.
#     "be" (data-driven, the eligible/primary k) marked. Anti-cherry-pick visual:
#     "be" sits on the plateau; fixed-k don't beat it.
# ===========================================================================
{
  # Build a k_sv-resolved table: prefer the winners table (has perm_typeI_mean,
  # loco_repro_scalar, ot_enrichment_or_up, k_sv, correction_id). Augment k_sv
  # from the manifest glob if winners lacks fixed-k rows.
  src <- NULL
  if (!is.null(W)) {
    keep <- intersect(c("cell_id","engine","correction_id","k_sv",
                        "perm_typeI_mean","loco_repro_scalar","ot_enrichment_or_up"),
                      names(W))
    src <- W[, keep, drop=FALSE]
  }
  # Manifest gives k_sv coverage even before metrics land (calibration/loco/ext
  # may still be NA -> their lines just won't draw).
  mf <- load_manifest()
  if (!is.null(mf)) {
    mf <- backfill_ids(mf)
    mfk <- mf[, intersect(c("cell_id","engine","correction_id","k_sv"), names(mf)), drop=FALSE]
    if (is.null(src)) src <- mfk
    else src <- dplyr::bind_rows(src, mfk[!mfk$cell_id %in% src$cell_id, , drop=FALSE])
  }
  cal <- load_calibration()
  if (!is.null(cal) && !is.null(src) && "perm_typeI_mean" %in% names(cal)) {
    cal2 <- cal[, intersect(c("cell_id","perm_typeI_mean"), names(cal)), drop=FALSE]
    if (!"perm_typeI_mean" %in% names(src)) src$perm_typeI_mean <- NA_real_
    src <- dplyr::rows_update(src, cal2, by="cell_id", unmatched="ignore")
  }

  metric_cols <- intersect(c("perm_typeI_mean","loco_repro_scalar","ot_enrichment_or_up"),
                           if (is.null(src)) character(0) else names(src))
  sv_present <- if (is.null(src)) FALSE else any(src$correction_id %in% SV_CORRECTIONS, na.rm=TRUE)
  has_kvar <- if (is.null(src)) FALSE else {
    kc <- src[src$correction_id %in% SV_CORRECTIONS, "k_sv", drop=TRUE]
    length(unique(kc[!is.na(kc) & kc != "na"])) >= 2
  }
  any_metric <- length(metric_cols) > 0 &&
    (if (is.null(src)) FALSE else any(sapply(metric_cols, function(c) any(is.finite(src[[c]])))))

  if (is.null(src) || !sv_present || !has_kvar) {
    skip("J4 sv sweep", "no SV/RUV correction with >=2 distinct k_sv values yet")
  } else if (!any_metric) {
    skip("J4 sv sweep", "k_sv grid present but no perm_typeI / loco / ot metric populated yet")
  } else {
    message("[render] J4 sv sweep")
    d4 <- src[src$correction_id %in% SV_CORRECTIONS, , drop=FALSE]
    # order k_sv on the canonical sweep axis
    klev <- c("be","2","5","10","20")
    d4 <- d4[d4$k_sv %in% klev, , drop=FALSE]
    d4$k_sv <- factor(d4$k_sv, levels=klev)
    metric_lab <- c(perm_typeI_mean="type-I error",
                    loco_repro_scalar="LOCO reproducibility",
                    ot_enrichment_or_up="OT enrichment OR (up)")
    long4 <- d4 %>%
      tidyr::pivot_longer(all_of(metric_cols), names_to="metric", values_to="value") %>%
      filter(is.finite(value)) %>%
      mutate(metric = factor(metric_lab[metric], levels=metric_lab[metric_cols]),
             corr   = pc(correction_id))
    # collapse engines: median per (correction, k_sv, metric)
    agg4 <- long4 %>%
      group_by(corr, correction_id, k_sv, metric) %>%
      summarise(value=median(value, na.rm=TRUE), .groups="drop")
    be_band <- function(df){ if ("type-I error" %in% levels(df$metric))
      data.frame(metric=factor("type-I error", levels=levels(df$metric))) else NULL }
    corr_pal <- setNames(
      colorRampPalette(c("#4aa2c2","#9b75d6","#e37faf"))(length(unique(agg4$corr))),
      sort(unique(agg4$corr)))
    pJ4 <- ggplot(agg4, aes(k_sv, value, group=corr, color=corr)) +
      geom_vline(xintercept=which(klev=="be"), linetype="dashed",
                 linewidth=0.3, color="#9E9E9E") +
      geom_line(linewidth=0.45) + geom_point(size=1.3) +
      scale_color_manual(values=corr_pal, name="correction") +
      facet_wrap(~metric, scales="free_y", nrow=1) +
      labs(title="Surrogate-variable sweep: does fixed-k beat data-driven 'be'?",
           subtitle="Median across engines vs number of surrogate variables; dashed = 'be' (data-driven, the eligible/primary k). Plateau at 'be' = no cherry-pick gain.",
           x="number of surrogate variables (k_sv)", y=NULL) +
      theme_masld()+theme_pub()+
      theme(legend.position="right", panel.spacing=unit(0.5,"lines"))
    save_fig(pJ4, file.path(OUT, "panelJ4_sv_sweep.pdf"),
             width=fig_full_width, height=3.2)
    mark(file.path(OUT, "panelJ4_sv_sweep.pdf"))
  }
}

# ===========================================================================
# J5. CONSENSUS CORE -- # genes called DEG (padj<0.05 & |logFC|>0.5) by exactly
#     k of the N ok cells, plus the >=k cumulative curve. Stable consensus core
#     vs method-unique tail. (panelI4 idea from figS_degx_multimethod.R.)
# ===========================================================================
{
  cell_files <- if (dir.exists(CELLS_DIR))
    list.files(CELLS_DIR, pattern="^cell_.*\\.csv$", full.names=TRUE) else character(0)
  if (length(cell_files) < 2) {
    skip("J5 consensus core", sprintf("need >=2 cells/cell_*.csv (found %d)", length(cell_files)))
  } else {
    message("[render] J5 consensus core (", length(cell_files), " cells)")
    deg_of <- function(f){
      d <- tryCatch(utils::read.csv(f, check.names=FALSE, stringsAsFactors=FALSE),
                    error=function(e) NULL)
      if (is.null(d) || !all(c("gene","logFC","padj") %in% names(d))) return(character(0))
      g <- d$gene[!is.na(d$padj) & d$padj < DEG_PADJ &
                  is.finite(d$logFC) & abs(d$logFC) > DEG_LFC]
      unique(as.character(g))
    }
    deg_sets <- lapply(cell_files, deg_of)
    N <- length(deg_sets)
    tab <- table(unlist(deg_sets))                # gene -> # cells calling it DEG
    if (length(tab) == 0) {
      skip("J5 consensus core", "no genes pass DEG = padj<0.05 & |logFC|>0.5 in any cell")
    } else {
      kcount <- as.integer(tab)
      exact <- as.data.frame(table(factor(kcount, levels=1:N)))
      names(exact) <- c("k","n_genes"); exact$k <- as.integer(as.character(exact$k))
      exact$atleast <- rev(cumsum(rev(exact$n_genes)))   # genes called by >= k cells
      lab_core <- sprintf("consensus core: %s genes in all %d cells",
                          format(exact$n_genes[exact$k==N], big.mark=","), N)
      lab_uniq <- sprintf("method-unique: %s genes in exactly 1",
                          format(exact$n_genes[exact$k==1], big.mark=","))
      bar_fill <- ifelse(exact$k==N, "#193c1e", ifelse(exact$k==1, "#9E9E9E", "#7fb069"))
      pJ5 <- ggplot(exact, aes(factor(k), n_genes)) +
        geom_col(aes(fill=I(bar_fill)), width=0.78, color="white", linewidth=0.25) +
        geom_line(aes(x=k, y=atleast), color="#C9265E", linewidth=0.5) +
        geom_point(aes(x=k, y=atleast), color="#C9265E", size=1.2) +
        annotate("text", x=N, y=max(exact$atleast), label=lab_core,
                 hjust=1, vjust=-0.6, size=PUB_GEOM_TEXT, color="grey20") +
        annotate("text", x=1, y=exact$n_genes[exact$k==1], label=lab_uniq,
                 hjust=0, vjust=-0.6, size=PUB_GEOM_TEXT, color="grey30") +
        scale_y_continuous(expand=expansion(mult=c(0,0.08)), labels=scales::comma) +
        labs(title="Consensus core vs method-unique tail",
             subtitle=sprintf("Bars = genes called DEG (padj<0.05 & |logFC|>0.5) by EXACTLY k of %d cells; magenta line = >= k (cumulative). Dark = all-cell core; grey = 1-cell unique.", N),
             x="number of cells calling the gene DEG", y="number of genes") +
        theme_masld()+theme_pub()+theme(legend.position="none")
      save_fig(pJ5, file.path(OUT, "panelJ5_consensus_core.pdf"),
               width=fig_col_width, height=3.4)
      mark(file.path(OUT, "panelJ5_consensus_core.pdf"))
    }
  }
}

# ---------------------------------------------------------------------------
message("\n[degx factorial] wrote ", length(written), " panel(s) to:\n  ", OUT)
invisible(lapply(written, function(p) cat("   -", basename(p), "\n")))
if (length(written) == 0)
  message("[degx factorial] no panels rendered -- all inputs absent/partial (expected on a fresh grid).")
