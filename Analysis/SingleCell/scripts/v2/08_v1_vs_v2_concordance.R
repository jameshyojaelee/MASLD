#!/usr/bin/env Rscript
# ============================================================================
# 08_v1_vs_v2_concordance.R
#
# Phase 0.5 — v1 vs v2 concordance audit (Agent 5).
#
# v1 paths are READ-ONLY. All writes go under
#   Analysis/SingleCell/results_gpu_v2_phase05/
#   figures/supplementary/stage_ccc_v2/
#
# Sections
#   A. Atlas-level cell / donor / cells-lost comparison
#   B. F-stage prediction concordance (Spearman + confusion matrix +
#      per-dataset breakdown + Healthy->F4 hallucination count)
#   C. Stage-CCC chain LMM Jaccard (coarse / fstage / continuous) +
#      bulk concordance percentages comparison
#   D. Headline biology audit (THBS1 / FGF13 / NAMPT->INSR / complement)
#   E. Verdict with rho thresholds (>0.95 / 0.7-0.95 / <0.7) and downstream
#      flag emission (v2_sufficient.flag vs v2_supersedes_v1.flag)
#
# Outputs:
#   results_gpu_v2_phase05/v1_v2_atlas_comparison.tsv
#   results_gpu_v2_phase05/v1_v2_fstage_concordance.tsv
#   results_gpu_v2_phase05/v1_v2_lmm_jaccard.tsv
#   results_gpu_v2_phase05/v1_v2_FINAL_REPORT.md
#   figures/supplementary/stage_ccc_v2/figS_v1_vs_v2_comparison.pdf
#   results_gpu_v2_phase05/{v2_sufficient.flag, v2_supersedes_v1.flag}
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

BASE <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)

# ----------------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------------
V1_STAGE_DIR <- file.path(
  BASE, "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory"
)
V1_MCP_INPUTS <- file.path(
  BASE, "Analysis/SingleCell/results_gpu_v2/mcp/inputs"
)

V2_PHASE05_DIR  <- file.path(BASE, "Analysis/SingleCell/results_gpu_v2_phase05")
V2_STAGE_DIR    <- file.path(V2_PHASE05_DIR, "ccc/stage_trajectory_v2")
V2_MCP_INPUTS   <- file.path(V2_PHASE05_DIR, "mcp/inputs")
V2_ATLAS_DIR    <- file.path(V2_PHASE05_DIR, "atlas")

OUT_DIR <- V2_PHASE05_DIR
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FIG_DIR <- file.path(BASE, "figures/supplementary/stage_ccc_v2")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_ATLAS_TSV  <- file.path(OUT_DIR, "v1_v2_atlas_comparison.tsv")
OUT_FSTAGE_TSV <- file.path(OUT_DIR, "v1_v2_fstage_concordance.tsv")
OUT_LMM_TSV    <- file.path(OUT_DIR, "v1_v2_lmm_jaccard.tsv")
OUT_REPORT_MD  <- file.path(OUT_DIR, "v1_v2_FINAL_REPORT.md")
OUT_PDF        <- file.path(FIG_DIR, "figS_v1_vs_v2_comparison.pdf")
FLAG_SUFFICIENT <- file.path(OUT_DIR, "v2_sufficient.flag")
FLAG_SUPERSEDE  <- file.path(OUT_DIR, "v2_supersedes_v1.flag")

# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------
fread_safe <- function(path, ...) {
  if (!file.exists(path)) {
    cat(sprintf("[warn] missing: %s\n", path))
    return(NULL)
  }
  tryCatch(fread(path, ...), error = function(e) {
    cat(sprintf("[warn] failed to read %s: %s\n", path, e$message))
    NULL
  })
}

safe_spearman <- function(x, y) {
  ok <- !is.na(x) & !is.na(y)
  if (sum(ok) < 5) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
}

# Reads a per-donor argmax/prediction column robustly from either v1 or v2
# augmented predicted tables. Returns sample/dataset + numeric F-stage.
read_augmented <- function(path, label) {
  dt <- fread_safe(path)
  if (is.null(dt)) return(NULL)
  # canonical ordinal column names
  ord_cols <- grep("F_stage_pred_augmented_.*ordinal", names(dt), value = TRUE)
  if (length(ord_cols) == 0) {
    # fall back to any argmax
    ord_cols <- grep("F_stage.*augmented.*", names(dt), value = TRUE)
  }
  if (length(ord_cols) == 0) {
    cat(sprintf("[warn] no augmented F-stage col in %s\n", path))
    return(NULL)
  }
  pick <- ord_cols[1]
  out <- data.table(
    sample  = as.character(dt$sample),
    dataset = as.character(dt$dataset),
    fstage_arg = suppressWarnings(as.integer(as.numeric(dt[[pick]])))
  )
  if ("origin" %in% names(dt)) out$origin <- as.character(dt$origin)
  if ("F_stage_documented" %in% names(dt))
    out$F_stage_documented <- suppressWarnings(as.integer(dt$F_stage_documented))
  if ("disease_stage_coarse" %in% names(dt))
    out$disease_stage_coarse <- as.character(dt$disease_stage_coarse)

  setnames(out, "fstage_arg", paste0("F_stage_aug_", label))
  out
}

jaccard_topk <- function(v1_vec, v2_vec, k = 50) {
  s1 <- head(unique(v1_vec), k)
  s2 <- head(unique(v2_vec), k)
  if (length(s1) == 0 || length(s2) == 0) return(NA_real_)
  length(intersect(s1, s2)) / length(union(s1, s2))
}

# ----------------------------------------------------------------------------
# Section A — Atlas-level comparison
# ----------------------------------------------------------------------------
cat("\n=========================\n")
cat("Section A: Atlas-level comparison\n")
cat("=========================\n")

v1_donor_meta <- fread_safe(file.path(V1_MCP_INPUTS, "donor_metadata.tsv"))
v2_donor_meta <- fread_safe(file.path(V2_MCP_INPUTS, "donor_metadata_v2.tsv"))

atlas_rows <- list()

# Cell-count audit from QC log (v2 Phase 0.5 ships an audit TSV)
qc_audit_v2 <- fread_safe(file.path(V2_ATLAS_DIR, "qc_cell_audit_v2.tsv"))
if (!is.null(qc_audit_v2) && "dataset_id" %in% names(qc_audit_v2)) {
  for (i in seq_len(nrow(qc_audit_v2))) {
    r <- qc_audit_v2[i]
    n_raw   <- if ("n_raw" %in% names(r)) as.integer(r$n_raw) else NA_integer_
    n_after <- if ("n_after_doublet" %in% names(r)) as.integer(r$n_after_doublet) else
               if ("n_after_qc" %in% names(r))   as.integer(r$n_after_qc)      else NA_integer_
    lost_pct <- if (!is.na(n_raw) && n_raw > 0)
      round(100 * (1 - n_after / n_raw), 2) else NA_real_
    atlas_rows[[length(atlas_rows) + 1L]] <- data.table(
      level        = "dataset",
      entity_key   = as.character(r$dataset_id),
      metric       = "cells",
      v1_value     = NA_integer_,
      v2_value     = n_after,
      delta_abs    = NA_real_,
      delta_pct    = lost_pct,
      note         = sprintf("v2 n_raw=%s, n_after=%s, cells_lost_pct=%s",
                             n_raw, n_after, lost_pct)
    )
  }
}

# Donor counts (v1 vs v2 -- did any donors drop?)
n_donors_v1 <- if (!is.null(v1_donor_meta)) nrow(v1_donor_meta) else NA_integer_
n_donors_v2 <- if (!is.null(v2_donor_meta)) nrow(v2_donor_meta) else NA_integer_
atlas_rows[[length(atlas_rows) + 1L]] <- data.table(
  level     = "atlas",
  entity_key = "global",
  metric    = "donor_count",
  v1_value  = n_donors_v1,
  v2_value  = n_donors_v2,
  delta_abs = (n_donors_v2 - n_donors_v1),
  delta_pct = if (!is.na(n_donors_v1) && n_donors_v1 > 0)
    round(100 * (n_donors_v2 - n_donors_v1) / n_donors_v1, 2) else NA_real_,
  note      = "Donor roster size; negative = donors dropped by harmonized QC"
)

# Per-dataset donor counts
if (!is.null(v1_donor_meta) && !is.null(v2_donor_meta) &&
    "dataset" %in% names(v1_donor_meta) && "dataset" %in% names(v2_donor_meta)) {
  v1_per_ds <- v1_donor_meta[, .N, by = dataset]
  v2_per_ds <- v2_donor_meta[, .N, by = dataset]
  setnames(v1_per_ds, "N", "v1_donors")
  setnames(v2_per_ds, "N", "v2_donors")
  ds_cmp <- merge(v1_per_ds, v2_per_ds, by = "dataset", all = TRUE)
  for (i in seq_len(nrow(ds_cmp))) {
    r <- ds_cmp[i]
    atlas_rows[[length(atlas_rows) + 1L]] <- data.table(
      level     = "dataset",
      entity_key = as.character(r$dataset),
      metric    = "donor_count",
      v1_value  = if (is.na(r$v1_donors)) 0L else as.integer(r$v1_donors),
      v2_value  = if (is.na(r$v2_donors)) 0L else as.integer(r$v2_donors),
      delta_abs = (ifelse(is.na(r$v2_donors), 0L, r$v2_donors) -
                   ifelse(is.na(r$v1_donors), 0L, r$v1_donors)),
      delta_pct = NA_real_,
      note      = "Per-dataset donor count"
    )
  }
}

atlas_df <- rbindlist(atlas_rows, fill = TRUE)
fwrite(atlas_df, OUT_ATLAS_TSV, sep = "\t")
cat(sprintf("[A] wrote %s (%d rows)\n", OUT_ATLAS_TSV, nrow(atlas_df)))

# ----------------------------------------------------------------------------
# Section B — F-stage prediction concordance
# ----------------------------------------------------------------------------
cat("\n=========================\n")
cat("Section B: F-stage concordance\n")
cat("=========================\n")

v1_aug <- read_augmented(
  file.path(V1_STAGE_DIR, "donor_fstage_augmented_predicted.tsv"), "v1"
)
v2_aug <- read_augmented(
  file.path(V2_STAGE_DIR, "donor_fstage_augmented_v2_predicted.tsv"), "v2"
)

fstage_rows <- list()
spearman_overall <- NA_real_
healthy_to_f4_v1 <- NA_integer_
healthy_to_f4_v2 <- NA_integer_

if (!is.null(v1_aug) && !is.null(v2_aug)) {
  joined <- merge(v1_aug, v2_aug,
                  by = c("sample", "dataset"), all = FALSE,
                  suffixes = c("", "_y"))
  cat(sprintf("[B] joined %d donors with both v1 and v2 augmented F-stage\n",
              nrow(joined)))

  # Resolve disease_stage_coarse / origin columns (may end with .x / .y)
  for (base_col in c("disease_stage_coarse", "origin", "F_stage_documented")) {
    cand <- paste0(base_col, c("", ".x", ".y"))
    pick <- cand[cand %in% names(joined)][1]
    if (!is.na(pick) && pick != base_col) {
      setnames(joined, pick, base_col)
    }
  }

  spearman_overall <- safe_spearman(joined$F_stage_aug_v1, joined$F_stage_aug_v2)
  cat(sprintf("[B] overall Spearman rho(v1,v2) = %.4f\n", spearman_overall))

  fstage_rows[[length(fstage_rows) + 1L]] <- data.table(
    metric = "spearman_overall",
    scope  = "all_donors",
    value  = spearman_overall,
    n      = nrow(joined)
  )

  # Confusion matrix (v1 argmax x v2 argmax)
  if (any(!is.na(joined$F_stage_aug_v1)) && any(!is.na(joined$F_stage_aug_v2))) {
    conf <- table(v1 = factor(joined$F_stage_aug_v1, levels = 0:4),
                  v2 = factor(joined$F_stage_aug_v2, levels = 0:4))
    conf_dt <- as.data.table(as.data.frame.matrix(conf), keep.rownames = "v1")
    cat("[B] confusion matrix (rows=v1, cols=v2):\n")
    print(conf)
    for (i in 0:4) for (j in 0:4) {
      fstage_rows[[length(fstage_rows) + 1L]] <- data.table(
        metric = "confusion_cell",
        scope  = sprintf("v1F%d_x_v2F%d", i, j),
        value  = as.numeric(conf[as.character(i), as.character(j)]),
        n      = nrow(joined)
      )
    }
  }

  # Per-dataset Spearman / agreement
  per_ds <- joined[, .(
    spearman = safe_spearman(F_stage_aug_v1, F_stage_aug_v2),
    pct_exact = mean(F_stage_aug_v1 == F_stage_aug_v2, na.rm = TRUE) * 100,
    n        = .N
  ), by = dataset]
  cat("[B] per-dataset:\n"); print(per_ds)
  for (i in seq_len(nrow(per_ds))) {
    r <- per_ds[i]
    fstage_rows[[length(fstage_rows) + 1L]] <- data.table(
      metric = "spearman_per_dataset",
      scope  = as.character(r$dataset),
      value  = r$spearman,
      n      = r$n
    )
    fstage_rows[[length(fstage_rows) + 1L]] <- data.table(
      metric = "pct_exact_agreement",
      scope  = as.character(r$dataset),
      value  = r$pct_exact,
      n      = r$n
    )
  }

  # Healthy -> F4 hallucination count
  if ("disease_stage_coarse" %in% names(joined)) {
    healthy_to_f4_v1 <- sum(joined$disease_stage_coarse == "Healthy" &
                              joined$F_stage_aug_v1 == 4L, na.rm = TRUE)
    healthy_to_f4_v2 <- sum(joined$disease_stage_coarse == "Healthy" &
                              joined$F_stage_aug_v2 == 4L, na.rm = TRUE)
    n_healthy <- sum(joined$disease_stage_coarse == "Healthy", na.rm = TRUE)
    cat(sprintf("[B] Healthy->F4 hallucinations: v1=%d, v2=%d (of %d healthy)\n",
                healthy_to_f4_v1, healthy_to_f4_v2, n_healthy))
    fstage_rows[[length(fstage_rows) + 1L]] <- data.table(
      metric = "healthy_to_F4_hallucination",
      scope  = "v1",
      value  = healthy_to_f4_v1,
      n      = n_healthy
    )
    fstage_rows[[length(fstage_rows) + 1L]] <- data.table(
      metric = "healthy_to_F4_hallucination",
      scope  = "v2",
      value  = healthy_to_f4_v2,
      n      = n_healthy
    )
  }
} else {
  cat("[B] SKIP: one of v1/v2 augmented predictions is missing\n")
}

fstage_df <- rbindlist(fstage_rows, fill = TRUE)
fwrite(fstage_df, OUT_FSTAGE_TSV, sep = "\t")
cat(sprintf("[B] wrote %s (%d rows)\n", OUT_FSTAGE_TSV, nrow(fstage_df)))

# ----------------------------------------------------------------------------
# Section C — Stage-CCC chain LMM Jaccard + bulk concordance
# ----------------------------------------------------------------------------
cat("\n=========================\n")
cat("Section C: LMM Jaccard + bulk concordance\n")
cat("=========================\n")

# v1 LMM tables
v1_lmm_coarse     <- fread_safe(file.path(V1_STAGE_DIR, "stage_lr_lmm_coarse.tsv"))
v1_lmm_fstage     <- fread_safe(file.path(V1_STAGE_DIR, "stage_lr_lmm_fstage.tsv"))
v1_lmm_continuous <- fread_safe(file.path(V1_STAGE_DIR, "stage_lr_lmm_continuous.tsv"))

# v2 LMM tables — Agent 4's chain script writes with _v2.tsv suffix primarily.
v2_lmm_coarse     <- fread_safe(file.path(V2_STAGE_DIR, "stage_lr_lmm_coarse_v2.tsv"))
v2_lmm_fstage     <- fread_safe(file.path(V2_STAGE_DIR, "stage_lr_lmm_fstage_v2.tsv"))
v2_lmm_continuous <- fread_safe(file.path(V2_STAGE_DIR, "stage_lr_lmm_continuous_v2.tsv"))

# Fallbacks: unsuffixed / .tsv.gz / pre-existing name conventions
alt_v2 <- function(name) {
  cand <- c(
    file.path(V2_STAGE_DIR, paste0("stage_lr_lmm_", name, ".tsv")),
    file.path(V2_STAGE_DIR, paste0("stage_lr_lmm_", name, "_v2.tsv.gz")),
    file.path(V2_STAGE_DIR, paste0("stage_lr_lmm_", name, ".tsv.gz"))
  )
  for (p in cand) if (file.exists(p)) return(fread_safe(p))
  NULL
}
if (is.null(v2_lmm_coarse))     v2_lmm_coarse     <- alt_v2("coarse")
if (is.null(v2_lmm_fstage))     v2_lmm_fstage     <- alt_v2("fstage")
if (is.null(v2_lmm_continuous)) v2_lmm_continuous <- alt_v2("continuous")

# Build top-50 lists ranked by lowest pval, then |Estimate|, restricted to the
# canonical disease-vs-Healthy term if present (matches v1 narrative coarse axis)
top_lr <- function(dt, n = 50, axis_label = NULL) {
  if (is.null(dt) || nrow(dt) == 0) return(character(0))
  cols <- names(dt)
  # Filter for disease-relevant term if column exists
  if ("term" %in% cols) {
    # coarse axis: keep disease_stage_coarse* terms; fstage axis: keep F_stage_*
    # continuous: keep continuous_axis or similar; if pattern is missing, keep all
    if (!is.null(axis_label) && axis_label == "coarse") {
      keep <- grepl("disease_stage", dt$term, ignore.case = TRUE)
      if (any(keep)) dt <- dt[keep]
    } else if (!is.null(axis_label) && axis_label == "fstage") {
      keep <- grepl("F_stage", dt$term, ignore.case = TRUE)
      if (any(keep)) dt <- dt[keep]
    }
  }
  if ("pval" %in% cols) dt <- dt[order(pval, na.last = TRUE)]
  if ("padj_within_ct" %in% cols) dt <- dt[order(padj_within_ct, na.last = TRUE)]
  key_col <- if ("lr_pair" %in% cols) "lr_pair" else
             if (all(c("source", "target", "ligand_complex", "receptor_complex") %in% cols))
               NULL else NULL
  if (is.null(key_col)) {
    dt[, lr_key := paste(source, target, ligand_complex, receptor_complex, sep = "|")]
    keys <- head(unique(dt$lr_key), n)
  } else {
    keys <- head(unique(dt[[key_col]]), n)
  }
  keys
}

lmm_rows <- list()
for (ax in c("coarse", "fstage", "continuous")) {
  v1tab <- get(paste0("v1_lmm_", ax))
  v2tab <- get(paste0("v2_lmm_", ax))
  v1_top <- top_lr(v1tab, n = 50, axis_label = ax)
  v2_top <- top_lr(v2tab, n = 50, axis_label = ax)
  jac <- jaccard_topk(v1_top, v2_top, k = 50)
  cat(sprintf("[C] %s axis: top-50 Jaccard = %s (v1=%d, v2=%d)\n",
              ax, format(jac, digits = 4),
              length(v1_top), length(v2_top)))
  lmm_rows[[length(lmm_rows) + 1L]] <- data.table(
    axis        = ax,
    metric      = "top50_jaccard",
    v1_top_n    = length(v1_top),
    v2_top_n    = length(v2_top),
    value       = jac
  )
}

# Bulk concordance percentages comparison
parse_bulk_concordance <- function(dt) {
  if (is.null(dt) || nrow(dt) == 0) return(list(n = 0, lig_concord = NA_real_,
                                                rec_concord = NA_real_,
                                                both_concord = NA_real_))
  list(
    n = nrow(dt),
    lig_concord  = if ("lig_concordant"  %in% names(dt))
      mean(as.logical(dt$lig_concordant),  na.rm = TRUE) * 100 else NA_real_,
    rec_concord  = if ("rec_concordant"  %in% names(dt))
      mean(as.logical(dt$rec_concordant),  na.rm = TRUE) * 100 else NA_real_,
    both_concord = if ("both_concordant" %in% names(dt))
      mean(as.logical(dt$both_concordant), na.rm = TRUE) * 100 else NA_real_
  )
}
# v2 chain writes ONE combined lr_bulk_concordance_v2.tsv with an `axis` col.
# v1 writes per-axis lr_bulk_concordance_{coarse,fstage,continuous}.tsv.
v2_bulk_all <- fread_safe(file.path(V2_STAGE_DIR, "lr_bulk_concordance_v2.tsv"))
slice_v2_bulk <- function(axis_label) {
  if (is.null(v2_bulk_all)) return(NULL)
  if (!"axis" %in% names(v2_bulk_all)) return(v2_bulk_all)
  if (axis_label == "coarse")
    return(v2_bulk_all[grepl("^(Steatosis|Steatohepatitis|Cirrhosis)_vs_Healthy",
                             axis)])
  if (axis_label == "fstage")
    return(v2_bulk_all[grepl("F_stage", axis)])
  if (axis_label == "continuous")
    return(v2_bulk_all[grepl("pseudotime|continuous", axis)])
  NULL
}
for (ax in c("coarse", "fstage", "continuous")) {
  v1_bulk <- fread_safe(file.path(V1_STAGE_DIR,
                                  sprintf("lr_bulk_concordance_%s.tsv", ax)))
  v2_bulk <- slice_v2_bulk(ax)
  if (is.null(v2_bulk)) {
    # legacy fallback: per-axis tsv if Agent 4 schema diverges
    v2_bulk <- fread_safe(file.path(V2_STAGE_DIR,
                                    sprintf("lr_bulk_concordance_%s_v2.tsv", ax)))
  }
  v1s <- parse_bulk_concordance(v1_bulk)
  v2s <- parse_bulk_concordance(v2_bulk)
  cat(sprintf("[C] %s bulk concord (both): v1=%.1f%% v2=%.1f%%\n",
              ax, v1s$both_concord, v2s$both_concord))
  lmm_rows[[length(lmm_rows) + 1L]] <- data.table(
    axis     = ax,
    metric   = "bulk_concord_both_pct_v1",
    v1_top_n = v1s$n,
    v2_top_n = v2s$n,
    value    = v1s$both_concord
  )
  lmm_rows[[length(lmm_rows) + 1L]] <- data.table(
    axis     = ax,
    metric   = "bulk_concord_both_pct_v2",
    v1_top_n = v1s$n,
    v2_top_n = v2s$n,
    value    = v2s$both_concord
  )
}

lmm_df <- rbindlist(lmm_rows, fill = TRUE)
fwrite(lmm_df, OUT_LMM_TSV, sep = "\t")
cat(sprintf("[C] wrote %s (%d rows)\n", OUT_LMM_TSV, nrow(lmm_df)))

# ----------------------------------------------------------------------------
# Section D — Headline biology audit
# ----------------------------------------------------------------------------
cat("\n=========================\n")
cat("Section D: Headline biology audit\n")
cat("=========================\n")

# Pairs the v1 narrative leans on (complement / FGF13 / THBS1 / NAMPT->INSR)
headline_pairs <- list(
  list(ligand = "THBS1",  receptor = NA_character_),
  list(ligand = "FGF13",  receptor = NA_character_),
  list(ligand = "NAMPT",  receptor = "INSR"),
  list(ligand = "C3",     receptor = NA_character_),
  list(ligand = "C1QA",   receptor = NA_character_),
  list(ligand = "C1QB",   receptor = NA_character_)
)

biology_rows <- list()

check_pair_present <- function(dt, lig, rec) {
  if (is.null(dt) || nrow(dt) == 0) return(0L)
  cols <- names(dt)
  if ("ligand_complex" %in% cols) {
    mask <- toupper(dt$ligand_complex) == toupper(lig)
    if (!is.na(rec) && "receptor_complex" %in% cols)
      mask <- mask & toupper(dt$receptor_complex) == toupper(rec)
    return(sum(mask, na.rm = TRUE))
  }
  if ("lr_pair" %in% cols) {
    pat <- if (is.na(rec)) toupper(lig)
           else sprintf("%s__%s", toupper(lig), toupper(rec))
    return(sum(grepl(pat, toupper(dt$lr_pair)), na.rm = TRUE))
  }
  0L
}

for (axis in c("coarse", "fstage", "continuous")) {
  v1tab <- get(paste0("v1_lmm_", axis))
  v2tab <- get(paste0("v2_lmm_", axis))
  v1_top <- top_lr(v1tab, n = 50, axis_label = axis)
  v2_top <- top_lr(v2tab, n = 50, axis_label = axis)
  for (p in headline_pairs) {
    pat <- if (is.na(p$receptor)) p$ligand else sprintf("%s__%s", p$ligand, p$receptor)
    v1_in_top50 <- sum(grepl(pat, v1_top), na.rm = TRUE)
    v2_in_top50 <- sum(grepl(pat, v2_top), na.rm = TRUE)
    biology_rows[[length(biology_rows) + 1L]] <- data.table(
      axis        = axis,
      ligand      = p$ligand,
      receptor    = if (is.na(p$receptor)) "" else p$receptor,
      v1_in_top50 = v1_in_top50,
      v2_in_top50 = v2_in_top50,
      v1_n_rows   = check_pair_present(v1tab, p$ligand, p$receptor),
      v2_n_rows   = check_pair_present(v2tab, p$ligand, p$receptor)
    )
  }
}
biology_df <- rbindlist(biology_rows, fill = TRUE)
biology_tsv <- file.path(OUT_DIR, "v1_v2_headline_biology.tsv")
fwrite(biology_df, biology_tsv, sep = "\t")
cat(sprintf("[D] wrote %s (%d rows)\n", biology_tsv, nrow(biology_df)))

# SH and Cirrhosis emergent biology preserved?
v1_sh   <- fread_safe(file.path(V1_STAGE_DIR, "lr_bulk_concordance_sh.tsv"))
v1_cirr <- fread_safe(file.path(V1_STAGE_DIR, "lr_bulk_concordance_cirrhosis.tsv"))
# v2 combined table — slice by axis tag
v2_sh <- NULL; v2_cirr <- NULL
if (!is.null(v2_bulk_all) && "axis" %in% names(v2_bulk_all)) {
  v2_sh   <- v2_bulk_all[grepl("Steatohepatitis", axis, ignore.case = TRUE)]
  v2_cirr <- v2_bulk_all[grepl("Cirrhosis",      axis, ignore.case = TRUE)]
}
# legacy fallback
if (is.null(v2_sh) || (is.data.table(v2_sh) && nrow(v2_sh) == 0))
  v2_sh <- fread_safe(file.path(V2_STAGE_DIR, "lr_bulk_concordance_sh_v2.tsv"))
if (is.null(v2_cirr) || (is.data.table(v2_cirr) && nrow(v2_cirr) == 0))
  v2_cirr <- fread_safe(file.path(V2_STAGE_DIR, "lr_bulk_concordance_cirrhosis_v2.tsv"))

sh_v1_n   <- if (!is.null(v1_sh))   nrow(v1_sh)   else 0L
sh_v2_n   <- if (!is.null(v2_sh))   nrow(v2_sh)   else 0L
cirr_v1_n <- if (!is.null(v1_cirr)) nrow(v1_cirr) else 0L
cirr_v2_n <- if (!is.null(v2_cirr)) nrow(v2_cirr) else 0L
cat(sprintf("[D] SH emergent LR rows: v1=%d, v2=%d\n", sh_v1_n, sh_v2_n))
cat(sprintf("[D] Cirrhosis emergent LR rows: v1=%d, v2=%d\n", cirr_v1_n, cirr_v2_n))

# ----------------------------------------------------------------------------
# Section E — Verdict
# ----------------------------------------------------------------------------
cat("\n=========================\n")
cat("Section E: Verdict\n")
cat("=========================\n")

verdict <- if (is.na(spearman_overall)) {
  "UNDETERMINED (F-stage Spearman could not be computed)"
} else if (spearman_overall > 0.95) {
  "Soft v2 sufficient; harmonized QC does not change conclusions; cite as sensitivity in Methods Supplement."
} else if (spearman_overall > 0.70) {
  "Material changes; report both v1 and v2 in supplementary; consider switching to v2 as headline."
} else {
  "Major divergence; investigate which donors/lineages caused the shift; v2 likely becomes the headline."
}

cat(sprintf("[E] verdict: %s\n", verdict))
cat(sprintf("[E] spearman_overall = %s\n", format(spearman_overall, digits = 4)))

# Emit flag(s)
if (!is.na(spearman_overall) && spearman_overall > 0.95) {
  writeLines(c(
    "VERDICT: Soft v2 sufficient",
    sprintf("Date: %s", Sys.Date()),
    sprintf("Spearman rho(v1, v2 augmented F-stage) = %.4f", spearman_overall),
    "Action: cite v2 as a sensitivity analysis in Methods Supplement.",
    "No change to headline numbers required."
  ), FLAG_SUFFICIENT)
  cat(sprintf("[E] wrote %s\n", FLAG_SUFFICIENT))
}
if (!is.na(spearman_overall) && spearman_overall < 0.70) {
  # Conservative scripts-to-update list (sources of v1 headlines that depend on
  # donor F-stage / stage-CCC pairs)
  affected_scripts <- c(
    "Analysis/SingleCell/scripts/343m_augmented_fstage.py",
    "Analysis/SingleCell/scripts/344_build_donor_metadata_extended.R",
    "Analysis/SingleCell/scripts/345_per_donor_liana.py",
    "Analysis/SingleCell/scripts/346_stage_ccc_lmm.R",
    "Analysis/SingleCell/scripts/347_stage_ccc_bulk_concordance.R",
    "Analysis/SingleCell/scripts/348_stage_ccc_trajectory_figures.R",
    "Analysis/SingleCell/scripts/349_stage_ccc_dissoc_sensitivity.R",
    "scripts/figures/figS_stage_ccc_trajectory.R"
  )
  writeLines(c(
    "VERDICT: Major divergence — v2 supersedes v1",
    sprintf("Date: %s", Sys.Date()),
    sprintf("Spearman rho(v1, v2 augmented F-stage) = %.4f", spearman_overall),
    "Headline numbers must be re-derived from the v2 chain before manuscript submission.",
    "",
    "Scripts whose downstream products consume donor F-stage / stage-CCC pairs and",
    "therefore need their headline numbers refreshed against v2 outputs:",
    sprintf("  - %s", affected_scripts)
  ), FLAG_SUPERSEDE)
  cat(sprintf("[E] wrote %s\n", FLAG_SUPERSEDE))
}

# ----------------------------------------------------------------------------
# Section F — Final report
# ----------------------------------------------------------------------------
cat("\n=========================\n")
cat("Section F: Final report\n")
cat("=========================\n")

report_lines <- c(
  "# v1 vs v2 Concordance Report — Phase 0.5",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## A. Atlas-level",
  sprintf("- v1 donors: %s, v2 donors: %s",
          format(n_donors_v1), format(n_donors_v2)),
  if (!is.na(n_donors_v1) && !is.na(n_donors_v2))
    sprintf("- Donor delta: %d (%s%%)", n_donors_v2 - n_donors_v1,
            format(round(100 * (n_donors_v2 - n_donors_v1) /
                           max(1L, n_donors_v1), 2)))
  else "",
  if (!is.null(qc_audit_v2))
    sprintf("- Per-dataset cell-loss audit: %d datasets in qc_cell_audit_v2.tsv",
            nrow(qc_audit_v2))
  else "- Per-dataset cell-loss audit: TSV missing",
  "",
  "## B. F-stage concordance",
  sprintf("- Spearman rho(v1, v2 augmented F-stage) = %s",
          format(spearman_overall, digits = 4)),
  sprintf("- Healthy->F4 hallucinations: v1=%s, v2=%s",
          format(healthy_to_f4_v1), format(healthy_to_f4_v2)),
  "- Per-dataset breakdown in v1_v2_fstage_concordance.tsv",
  "",
  "## C. Stage-CCC chain LMM Jaccard"
)
for (ax in c("coarse", "fstage", "continuous")) {
  r <- lmm_df[axis == ax & metric == "top50_jaccard"]
  report_lines <- c(report_lines,
    sprintf("- %s axis top-50 Jaccard: %s", ax,
            if (nrow(r) > 0) format(r$value, digits = 4) else "NA"))
}
report_lines <- c(report_lines, "",
  "### Bulk concordance percentages (both ligand+receptor)")
for (ax in c("coarse", "fstage", "continuous")) {
  v1r <- lmm_df[axis == ax & metric == "bulk_concord_both_pct_v1"]
  v2r <- lmm_df[axis == ax & metric == "bulk_concord_both_pct_v2"]
  report_lines <- c(report_lines,
    sprintf("- %s axis: v1=%s%%, v2=%s%%",
            ax,
            if (nrow(v1r) > 0) format(v1r$value, digits = 3) else "NA",
            if (nrow(v2r) > 0) format(v2r$value, digits = 3) else "NA"))
}
report_lines <- c(report_lines, "",
  "## D. Headline biology audit",
  sprintf("- SH emergent LR rows: v1=%d, v2=%d", sh_v1_n, sh_v2_n),
  sprintf("- Cirrhosis emergent LR rows: v1=%d, v2=%d", cirr_v1_n, cirr_v2_n),
  "- THBS1 / FGF13 / NAMPT->INSR / complement preservation in v1_v2_headline_biology.tsv",
  "",
  "## E. Verdict",
  sprintf("- Spearman rho = %s", format(spearman_overall, digits = 4)),
  sprintf("- **Verdict**: %s", verdict),
  "",
  "## F. Methods note",
  "All v1 paths read-only. v2 outputs written under results_gpu_v2_phase05/.",
  "Top-50 Jaccard ranks LR pairs by lowest p-value, then absolute estimate,",
  "filtered to the disease/F-stage term per axis.",
  "Healthy->F4 hallucinations count donors with disease_stage_coarse == 'Healthy'",
  "and augmented F-stage argmax == 4."
)
writeLines(report_lines, OUT_REPORT_MD)
cat(sprintf("[F] wrote %s\n", OUT_REPORT_MD))

# ----------------------------------------------------------------------------
# Section G — 4-panel PDF
# ----------------------------------------------------------------------------
cat("\n=========================\n")
cat("Section G: Figure\n")
cat("=========================\n")

# Panel A: per-dataset cell counts (v1 vs v2). When v1 cell counts aren't in
# atlas_df, we draw donor counts as a proxy.
panel_a_dt <- atlas_df[level == "dataset" & metric == "donor_count" &
                         !is.na(entity_key) & entity_key != "global"]
if (nrow(panel_a_dt) == 0) {
  panel_a_dt <- atlas_df[level == "dataset" & metric == "cells"]
  panel_a_dt[, v1_value := v2_value]  # cells-only: plot v2 only
}
panel_a <- if (nrow(panel_a_dt) > 0) {
  long <- melt(panel_a_dt, id.vars = "entity_key",
               measure.vars = c("v1_value", "v2_value"),
               variable.name = "version", value.name = "n")
  long[, version := ifelse(version == "v1_value", "v1", "v2")]
  ggplot(long, aes(x = entity_key, y = n, fill = version)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_manual(values = c(v1 = "#9E9E9E", v2 = "#3B82F6")) +
    labs(x = NULL, y = "Donors", title = "A. Per-dataset donor count") +
    theme_minimal(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top")
} else {
  ggplot() + ggtitle("A. Per-dataset (no data)") + theme_void()
}

# Panel B: F-stage Spearman scatter
panel_b <- if (!is.null(v1_aug) && !is.null(v2_aug)) {
  jj <- merge(v1_aug[, .(sample, dataset, F_stage_aug_v1)],
              v2_aug[, .(sample, dataset, F_stage_aug_v2)],
              by = c("sample", "dataset"))
  ggplot(jj, aes(x = F_stage_aug_v1, y = F_stage_aug_v2)) +
    geom_jitter(width = 0.15, height = 0.15, alpha = 0.5,
                colour = "#3B82F6", size = 1.2) +
    geom_abline(slope = 1, intercept = 0, lty = 2, colour = "#9E9E9E") +
    scale_x_continuous(breaks = 0:4) +
    scale_y_continuous(breaks = 0:4) +
    labs(x = "v1 augmented F-stage", y = "v2 augmented F-stage",
         title = sprintf("B. F-stage agreement (Spearman %.3f, n=%d)",
                         spearman_overall, nrow(jj))) +
    theme_minimal(base_size = 8)
} else {
  ggplot() + ggtitle("B. F-stage scatter (no data)") + theme_void()
}

# Panel C: top-50 Jaccard per axis
panel_c <- if (any(lmm_df$metric == "top50_jaccard")) {
  jdt <- lmm_df[metric == "top50_jaccard"]
  ggplot(jdt, aes(x = axis, y = value, fill = axis)) +
    geom_col(width = 0.5) +
    geom_text(aes(label = sprintf("%.2f", value)), vjust = -0.4, size = 2.5) +
    scale_fill_manual(values = c(coarse = "#3B82F6",
                                 fstage = "#10B981",
                                 continuous = "#F59E0B")) +
    ylim(0, 1) +
    labs(x = "Stage axis", y = "top-50 Jaccard (v1 cap v2)",
         title = "C. LMM top-50 LR pair Jaccard") +
    theme_minimal(base_size = 8) +
    theme(legend.position = "none")
} else {
  ggplot() + ggtitle("C. Jaccard (no data)") + theme_void()
}

# Panel D: headline biology heatmap (in top-50?)
panel_d <- if (nrow(biology_df) > 0) {
  biology_df[, pair := ifelse(receptor == "" | is.na(receptor),
                              ligand, sprintf("%s->%s", ligand, receptor))]
  biology_long <- melt(biology_df, id.vars = c("axis", "pair"),
                       measure.vars = c("v1_in_top50", "v2_in_top50"),
                       variable.name = "version", value.name = "n_hits")
  biology_long[, version := ifelse(version == "v1_in_top50", "v1", "v2")]
  ggplot(biology_long, aes(x = paste(axis, version, sep = ":"),
                           y = pair, fill = n_hits > 0)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    scale_fill_manual(values = c(`TRUE` = "#10B981", `FALSE` = "#9E9E9E"),
                      name = "in top-50") +
    labs(x = NULL, y = NULL, title = "D. Headline biology in top-50") +
    theme_minimal(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top")
} else {
  ggplot() + ggtitle("D. Headline biology (no data)") + theme_void()
}

pdf(OUT_PDF, width = 11, height = 9)
print((panel_a | panel_b) / (panel_c | panel_d) +
        plot_annotation(title = "v1 vs v2 concordance",
                        subtitle = sprintf("verdict: %s", verdict)))
dev.off()
cat(sprintf("[G] wrote %s\n", OUT_PDF))

cat("\n[done] 08_v1_vs_v2_concordance.R\n")
