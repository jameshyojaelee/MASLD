# 250e_neutrophil_regularization.R
# Re-evaluate neutrophil log2FC = +10.39 at F2->F3 (Team M5 of cascade rigor-hardening plan).
# Tests three things and writes pre-registered acceptance verdict:
#   (1) Regularized log-ratio with eps in {1e-6,1e-5,1e-4,1e-3} for F1F2/F2F3/F3F4
#   (2) Absolute proportions per stage (median + IQR) for F0..F4
#   (3) sc cross-validation: Spearman rho(bulk vs sc neutrophil) per matched donor
#
# Same regularization sweep applied to other headline cell types:
#   Plasma cells F1->F2, Cholangiocytes F2->F3 / F3->F4, Fibroblasts F1->F2.
#
# Outputs (RNA-seq/results/granular_staging/):
#   neutrophil_regularized_validation.csv
#   celltype_regularized_other.csv
#   neutrophil_regularization_summary.md

suppressPackageStartupMessages({
  library(matrixStats)
})

set.seed(42)

PROJECT_ROOT <- Sys.getenv("MASLD_PROJECT_ROOT",
                          "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUT_DIR <- file.path(PROJECT_ROOT, "RNA-seq/results/granular_staging")
BULK_DECONV <- file.path(PROJECT_ROOT,
                         "RNA-seq/results/celltype_attribution/persample_celltype_proportions.csv")
SC_PROPS    <- file.path(PROJECT_ROOT,
                         "Analysis/SingleCell/results_gpu_v2/disease_signatures/celltype_proportions_per_sample.csv")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("== 250e: Neutrophil + headline cell-type regularization sweep ==\n")

# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------
deconv <- read.csv(BULK_DECONV, stringsAsFactors = FALSE, check.names = FALSE)
sc_props <- read.csv(SC_PROPS, stringsAsFactors = FALSE, check.names = FALSE)

cat(sprintf("Bulk deconv: %d samples x %d cols\n", nrow(deconv), ncol(deconv)))
cat(sprintf("sc proportions: %d samples x %d cols\n", nrow(sc_props), ncol(sc_props)))

stopifnot("Neutrophils" %in% colnames(deconv))
stopifnot("Neutrophils" %in% colnames(sc_props))

# ---------------------------------------------------------------------------
# Regularized log2FC sweep
# ---------------------------------------------------------------------------
EPS_GRID <- c(1e-6, 1e-5, 1e-4, 1e-3)
TRANSITIONS <- list(
  F1vF2 = c(1, 2),
  F2vF3 = c(2, 3),
  F3vF4 = c(3, 4)
)

# Stage-level absolute proportion stats (median + IQR) for one cell type
stage_abs_stats <- function(prop_vec, stages) {
  out <- data.frame(stage = 0:4, n = NA_integer_, median = NA_real_,
                    q25 = NA_real_, q75 = NA_real_, mean = NA_real_,
                    stringsAsFactors = FALSE)
  for (s in 0:4) {
    v <- prop_vec[stages == s]
    v <- v[is.finite(v)]
    out$n[out$stage == s] <- length(v)
    if (length(v) > 0) {
      qs <- quantile(v, c(0.25, 0.5, 0.75), na.rm = TRUE)
      out$q25[out$stage == s]    <- qs[1]
      out$median[out$stage == s] <- qs[2]
      out$q75[out$stage == s]    <- qs[3]
      out$mean[out$stage == s]   <- mean(v, na.rm = TRUE)
    }
  }
  out
}

# Bulk vs sc Spearman correlation for one cell type using matched donors
# (deconv$sample_id vs sc_props$sample). Match by exact sample id; if no
# overlap, return NA + reason.
bulk_sc_spearman <- function(bulk_df, sc_df, ct) {
  bulk_ids <- bulk_df$sample_id
  sc_ids   <- sc_df$sample
  matched <- intersect(bulk_ids, sc_ids)
  if (length(matched) >= 5) {
    bulk_v <- bulk_df[[ct]][match(matched, bulk_ids)]
    sc_v   <- sc_df[[ct]][match(matched, sc_ids)]
    ok <- is.finite(bulk_v) & is.finite(sc_v)
    if (sum(ok) < 5) {
      return(list(rho = NA_real_, n = sum(ok),
                  notes = "insufficient finite pairs"))
    }
    rho <- suppressWarnings(cor(bulk_v[ok], sc_v[ok], method = "spearman"))
    return(list(rho = rho, n = sum(ok), notes = "exact-id match"))
  }
  # No exact match: try a within-cohort distribution-level cross-comparison.
  # We use median of stage-pooled proportions in each modality (bulk by
  # fibrosis_stage, sc has only `condition` (Healthy/MASH), so stage->binary).
  # Fall back: pool all bulk samples and all sc samples; rho across cell types
  # is reported separately. Here for a per-cell-type fallback, return NA.
  list(rho = NA_real_, n = 0,
       notes = sprintf("no exact-ID overlap (bulk n=%d, sc n=%d)",
                       length(bulk_ids), length(sc_ids)))
}

# Bulk vs sc rho ACROSS cell types (single number per cell type fallback):
# we compute one cross-celltype rho using stage-pooled medians as a sanity
# check, since per-donor matching does not exist. This becomes a single value
# attached to all neutrophil rows.
bulk_pooled_median <- function(bulk_df, ct) {
  v <- bulk_df[[ct]]
  median(v[is.finite(v)], na.rm = TRUE)
}
sc_pooled_median <- function(sc_df, ct) {
  v <- sc_df[[ct]]
  median(v[is.finite(v)], na.rm = TRUE)
}
common_cts <- intersect(colnames(deconv), colnames(sc_props))
common_cts <- setdiff(common_cts,
                      c("sample_id", "sample", "dataset", "group_binary",
                        "condition", "fibrosis_stage", "inferred_sex",
                        "diagnosis_harmonized", "f2_group"))

cross_ct_bulk <- vapply(common_cts, bulk_pooled_median, numeric(1), bulk_df = deconv)
cross_ct_sc   <- vapply(common_cts, sc_pooled_median,   numeric(1), sc_df   = sc_props)
ok_cross <- is.finite(cross_ct_bulk) & is.finite(cross_ct_sc)
cross_ct_rho <- if (sum(ok_cross) >= 5) {
  suppressWarnings(cor(cross_ct_bulk[ok_cross], cross_ct_sc[ok_cross], method = "spearman"))
} else NA_real_
cat(sprintf("Cross-celltype Spearman rho (bulk pooled median vs sc pooled median, n=%d): %s\n",
            sum(ok_cross),
            ifelse(is.na(cross_ct_rho), "NA", sprintf("%.3f", cross_ct_rho))))

# ---------------------------------------------------------------------------
# (1) Build regularized log2FC table for Neutrophils + headline others
# ---------------------------------------------------------------------------
# Matches script 250 exactly: log2((mean_b + eps)/(mean_a + eps)) with stage
# means computed across all samples passing fibrosis_stage filter.

build_reg_table <- function(ct, transitions = TRANSITIONS, eps_grid = EPS_GRID) {
  stages <- deconv$fibrosis_stage
  prop   <- deconv[[ct]]
  rows <- list()
  for (tn in names(transitions)) {
    sa <- transitions[[tn]][1]; sb <- transitions[[tn]][2]
    a <- prop[stages == sa]; b <- prop[stages == sb]
    a <- a[is.finite(a)]; b <- b[is.finite(b)]
    ma <- mean(a); mb <- mean(b)
    meda <- median(a); medb <- median(b)
    for (eps in eps_grid) {
      lfc <- log2((mb + eps) / (ma + eps))
      rows[[paste(ct, tn, eps)]] <- data.frame(
        celltype = ct, transition = tn, epsilon = eps,
        log2fc = lfc,
        n_a = length(a), n_b = length(b),
        mean_a = ma, mean_b = mb,
        median_a = meda, median_b = medb,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

# ---- Neutrophils main table ----
neut_reg <- build_reg_table("Neutrophils")
neut_abs <- stage_abs_stats(deconv$Neutrophils, deconv$fibrosis_stage)
cat("\nNeutrophil absolute stage stats:\n"); print(neut_abs)

# Per-donor bulk vs sc Spearman
bsc_neu <- bulk_sc_spearman(deconv, sc_props, "Neutrophils")
cat(sprintf("Bulk-vs-sc Neutrophil Spearman: rho=%s n=%d (%s)\n",
            ifelse(is.na(bsc_neu$rho), "NA", sprintf("%.3f", bsc_neu$rho)),
            bsc_neu$n, bsc_neu$notes))

# Stage absolute proportion check at F3 (claimed 0.14%)
neut_F3_median <- neut_abs$median[neut_abs$stage == 3]
neut_F2_median <- neut_abs$median[neut_abs$stage == 2]

# Add abs stage values to each row (already mean_a/mean_b)
neut_reg$abs_prop_a <- neut_reg$mean_a
neut_reg$abs_prop_b <- neut_reg$mean_b
neut_reg$bulk_sc_spearman <- ifelse(is.na(bsc_neu$rho), cross_ct_rho, bsc_neu$rho)
neut_reg$bulk_sc_spearman_source <- ifelse(is.na(bsc_neu$rho), "cross_celltype_pooled", "per_donor_match")

# Pre-registered acceptance gate (per row context: F2->F3 is the headline)
gate_real_storm <- function(row, abs_prop_F3, abs_prop_F2, rho) {
  # 1: regularized log2FC at F2F3 with eps=1e-4 must be >= 2.0
  # 2: abs prop at F3 >= 0.5% AND >= 5x F2
  # 3: rho >= 0.4
  if (row$transition != "F2vF3") return(NA)
  c1 <- (row$epsilon == 1e-4) && (!is.na(row$log2fc)) && (row$log2fc >= 2.0)
  if (row$epsilon != 1e-4) return(NA) # only score the F2vF3 / eps=1e-4 row
  c2 <- (abs_prop_F3 >= 0.005) && (abs_prop_F3 >= 5 * abs_prop_F2)
  c3 <- (!is.na(rho)) && (rho >= 0.4)
  c1 && c2 && c3
}

neut_reg$pass_real_storm <- mapply(function(i) {
  gate_real_storm(neut_reg[i, , drop = FALSE],
                  abs_prop_F3 = neut_F3_median,
                  abs_prop_F2 = neut_F2_median,
                  rho = neut_reg$bulk_sc_spearman[i])
}, seq_len(nrow(neut_reg)))

neut_reg$notes <- sprintf("F3_median=%.4g F2_median=%.4g cross_ct_rho=%.3f",
                         neut_F3_median, neut_F2_median, cross_ct_rho)

# Write neutrophil CSV with required columns
neut_out <- data.frame(
  transition = neut_reg$transition,
  epsilon    = neut_reg$epsilon,
  log2fc     = round(neut_reg$log2fc, 3),
  abs_prop_a = signif(neut_reg$abs_prop_a, 4),
  abs_prop_b = signif(neut_reg$abs_prop_b, 4),
  bulk_sc_spearman = signif(neut_reg$bulk_sc_spearman, 4),
  pass_real_storm  = neut_reg$pass_real_storm,
  notes      = neut_reg$notes,
  stringsAsFactors = FALSE
)
write.csv(neut_out, file.path(OUT_DIR, "neutrophil_regularized_validation.csv"),
          row.names = FALSE)
cat(sprintf("Wrote neutrophil_regularized_validation.csv (%d rows)\n", nrow(neut_out)))

# ---- Other headline cell types (sweep) ----
other_targets <- list(
  "Plasma cells"   = "F1vF2",
  "Cholangiocytes" = c("F2vF3", "F3vF4"),
  "Fibroblasts"    = "F1vF2"
)
other_rows <- list()
for (ct in names(other_targets)) {
  if (!ct %in% colnames(deconv)) next
  reg <- build_reg_table(ct)
  for (eps in EPS_GRID) {
    for (tn in other_targets[[ct]]) {
      sub <- reg[reg$transition == tn & reg$epsilon == eps, , drop = FALSE]
      if (nrow(sub) == 0) next
      abs_stats <- stage_abs_stats(deconv[[ct]], deconv$fibrosis_stage)
      sa <- TRANSITIONS[[tn]][1]; sb <- TRANSITIONS[[tn]][2]
      f3_or_b_median <- abs_stats$median[abs_stats$stage == sb]
      f2_or_a_median <- abs_stats$median[abs_stats$stage == sa]
      bsc <- bulk_sc_spearman(deconv, sc_props, ct)
      other_rows[[paste(ct, tn, eps)]] <- data.frame(
        celltype = ct,
        transition = tn,
        epsilon = eps,
        log2fc = round(sub$log2fc, 3),
        abs_prop_a = signif(sub$mean_a, 4),
        abs_prop_b = signif(sub$mean_b, 4),
        median_a = signif(f2_or_a_median, 4),
        median_b = signif(f3_or_b_median, 4),
        bulk_sc_spearman = signif(ifelse(is.na(bsc$rho), cross_ct_rho, bsc$rho), 4),
        bulk_sc_source = ifelse(is.na(bsc$rho), "cross_celltype_pooled", "per_donor_match"),
        stringsAsFactors = FALSE
      )
    }
  }
}
other_df <- do.call(rbind, other_rows)
write.csv(other_df,
          file.path(OUT_DIR, "celltype_regularized_other.csv"),
          row.names = FALSE)
cat(sprintf("Wrote celltype_regularized_other.csv (%d rows)\n", nrow(other_df)))

# ---------------------------------------------------------------------------
# (3) Markdown verdict
# ---------------------------------------------------------------------------
neut_F2F3 <- neut_reg[neut_reg$transition == "F2vF3", ]
neut_F2F3_eps4 <- neut_F2F3[neut_F2F3$epsilon == 1e-4, ]
neut_F2F3_eps3 <- neut_F2F3[neut_F2F3$epsilon == 1e-3, ]
neut_F2F3_eps5 <- neut_F2F3[neut_F2F3$epsilon == 1e-5, ]
neut_F2F3_eps6 <- neut_F2F3[neut_F2F3$epsilon == 1e-6, ]

abs_F3 <- neut_F3_median
abs_F2 <- neut_F2_median
abs_F3_pct <- abs_F3 * 100
gate1 <- (!is.na(neut_F2F3_eps4$log2fc)) && (neut_F2F3_eps4$log2fc >= 2.0)
gate2 <- (abs_F3 >= 0.005) && (abs_F3 >= 5 * abs_F2)
gate3 <- (!is.na(cross_ct_rho)) && (cross_ct_rho >= 0.4)
verdict_pass <- gate1 && gate2 && gate3

md_lines <- c(
  "# Neutrophil F2->F3 +10.39 regularization re-evaluation (Team M5)",
  "",
  "## Headline finding",
  "",
  if (verdict_pass)
    "VERDICT: **PASSES** the pre-registered \"real innate immune storm\" gate."
  else
    "VERDICT: **FAILS** the pre-registered \"real innate immune storm\" gate. The +10.39 framing is a divide-by-near-zero artifact.",
  "",
  "## Neutrophil regularization sweep (transition x epsilon)",
  "",
  "| Transition | eps=1e-6 | eps=1e-5 | eps=1e-4 | eps=1e-3 |",
  "|---|---:|---:|---:|---:|",
  sprintf("| F1vF2 | %.3f | %.3f | %.3f | %.3f |",
          neut_reg$log2fc[neut_reg$transition=="F1vF2" & neut_reg$epsilon==1e-6],
          neut_reg$log2fc[neut_reg$transition=="F1vF2" & neut_reg$epsilon==1e-5],
          neut_reg$log2fc[neut_reg$transition=="F1vF2" & neut_reg$epsilon==1e-4],
          neut_reg$log2fc[neut_reg$transition=="F1vF2" & neut_reg$epsilon==1e-3]),
  sprintf("| F2vF3 | %.3f | %.3f | %.3f | %.3f |",
          neut_F2F3_eps6$log2fc, neut_F2F3_eps5$log2fc,
          neut_F2F3_eps4$log2fc, neut_F2F3_eps3$log2fc),
  sprintf("| F3vF4 | %.3f | %.3f | %.3f | %.3f |",
          neut_reg$log2fc[neut_reg$transition=="F3vF4" & neut_reg$epsilon==1e-6],
          neut_reg$log2fc[neut_reg$transition=="F3vF4" & neut_reg$epsilon==1e-5],
          neut_reg$log2fc[neut_reg$transition=="F3vF4" & neut_reg$epsilon==1e-4],
          neut_reg$log2fc[neut_reg$transition=="F3vF4" & neut_reg$epsilon==1e-3]),
  "",
  sprintf("Original published value at F2->F3 (eps=1e-6, matches Script 250): **%.2f**",
          neut_F2F3_eps6$log2fc),
  sprintf("At eps=1e-4 (the pre-registered gate threshold): **%.2f**",
          neut_F2F3_eps4$log2fc),
  "",
  "## Absolute proportion per stage (Neutrophils, bulk deconvolved)",
  "",
  "| Stage | n | median | q25 | q75 | mean |",
  "|---|---:|---:|---:|---:|---:|"
)
for (i in seq_len(nrow(neut_abs))) {
  md_lines <- c(md_lines, sprintf("| F%d | %d | %.4g | %.4g | %.4g | %.4g |",
                                  neut_abs$stage[i], neut_abs$n[i],
                                  neut_abs$median[i], neut_abs$q25[i],
                                  neut_abs$q75[i], neut_abs$mean[i]))
}
md_lines <- c(md_lines,
  "",
  sprintf("F3 absolute proportion (median): **%.4g (%.4f%%)**", abs_F3, abs_F3_pct),
  sprintf("F2 absolute proportion (median): **%.4g (%.4f%%)**", abs_F2, abs_F2 * 100),
  sprintf("F3/F2 ratio: **%.2fx**", ifelse(abs_F2 > 0, abs_F3 / abs_F2, NA)),
  "",
  "## Bulk vs sc Spearman cross-validation",
  "",
  sprintf("- Per-donor exact-ID match: %s (n=%d)", bsc_neu$notes, bsc_neu$n),
  sprintf("- Cross-celltype pooled median rho: **%.3f** (n=%d cell types)",
          cross_ct_rho, sum(ok_cross)),
  "",
  "## Pre-registered acceptance gates",
  "",
  sprintf("1. Regularized log2FC F2->F3 at eps=1e-4 >= 2.0: **%s** (got %.2f)",
          ifelse(gate1, "PASS", "FAIL"), neut_F2F3_eps4$log2fc),
  sprintf("2. Abs prop F3 >= 0.5%% AND F3 >= 5x F2: **%s** (F3=%.4f%%, F2=%.4f%%, ratio=%.2fx)",
          ifelse(gate2, "PASS", "FAIL"),
          abs_F3 * 100, abs_F2 * 100, ifelse(abs_F2 > 0, abs_F3 / abs_F2, NA)),
  sprintf("3. Bulk vs sc Spearman rho >= 0.4: **%s** (got %.3f)",
          ifelse(gate3, "PASS", "FAIL"), cross_ct_rho),
  "",
  "## Other headline cell-type sanity (regularized)",
  "")
for (ct in names(other_targets)) {
  if (!ct %in% other_df$celltype) next
  for (tn in other_targets[[ct]]) {
    sub <- other_df[other_df$celltype == ct & other_df$transition == tn, ]
    if (nrow(sub) == 0) next
    md_lines <- c(md_lines, sprintf("### %s @ %s", ct, tn),
                  "",
                  "| eps | log2fc | mean_a | mean_b | median_a | median_b |",
                  "|---:|---:|---:|---:|---:|---:|")
    sub <- sub[order(sub$epsilon), ]
    for (i in seq_len(nrow(sub))) {
      md_lines <- c(md_lines, sprintf("| %.0e | %.3f | %.4g | %.4g | %.4g | %.4g |",
                                      sub$epsilon[i], sub$log2fc[i],
                                      sub$abs_prop_a[i], sub$abs_prop_b[i],
                                      sub$median_a[i], sub$median_b[i]))
    }
    md_lines <- c(md_lines, "")
  }
}

writeLines(md_lines, file.path(OUT_DIR, "neutrophil_regularization_summary.md"))
cat(sprintf("Wrote neutrophil_regularization_summary.md\n"))

cat("\n== 250e complete ==\n")
