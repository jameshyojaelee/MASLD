#!/usr/bin/env Rscript
# ============================================================================
# 07_chain_346_349_v2.R
#
# Phase 0.5 v2 equivalent of Scripts 346 / 347 / 348 / 349.
#
# Runs in one driver because the four steps share a per-donor LIANA long-form
# input table and a donor_metadata_v2.tsv covariate joins.
#
# Inputs (all read-only on v1; v2 is read-write):
#   - results_gpu_v2_phase05/ccc/stage_trajectory_v2/per_donor_lr_v2/*.parquet
#   - results_gpu_v2_phase05/mcp/inputs/donor_metadata_v2.tsv
#   - RNA-seq/Human/.../dream_results_ashr.csv         (bulk; read-only)
#   - RNA-seq/Human/.../dream_results_stage_*.csv      (bulk stage anchors; read-only)
#
# Outputs (all under results_gpu_v2_phase05):
#   - stage_lr_lmm_coarse_v2.tsv
#   - stage_lr_lmm_fstage_v2.tsv
#   - stage_lr_lmm_continuous_v2.tsv
#   - loo_replication_rate_per_lr_v2.tsv
#   - lr_bulk_concordance_v2.tsv
#   - figures/supplementary/stage_ccc_v2/figS_stage_ccc_trajectory_v2_atlas.pdf
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
  library(lmerTest)
  library(parallel)
  library(ggplot2)
  library(patchwork)
})

BASE <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)

N_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "16"))
options(mc.cores = N_CORES)
cat(sprintf("[parallel] using %d cores\n", N_CORES))

V2_STAGE_DIR <- file.path(
  BASE, "Analysis/SingleCell/results_gpu_v2_phase05/ccc/stage_trajectory_v2"
)
PER_DONOR_DIR <- file.path(V2_STAGE_DIR, "per_donor_lr_v2")
META_V2 <- file.path(
  BASE,
  "Analysis/SingleCell/results_gpu_v2_phase05/mcp/inputs/donor_metadata_v2.tsv"
)
MERGED_TSV <- file.path(V2_STAGE_DIR, "all_donor_lr_scores_v2.tsv.gz")
FIG_DIR <- file.path(BASE, "figures/supplementary/stage_ccc_v2")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(V2_STAGE_DIR, recursive = TRUE, showWarnings = FALSE)

INT_DIR <- file.path(
  BASE,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration"
)
BULK_ALL <- file.path(INT_DIR, "dream_results_ashr.csv")
BULK_STEATOSIS <- file.path(INT_DIR, "dream_results_stage_steatosis.csv")
BULK_SH <- file.path(INT_DIR, "dream_results_stage_sh.csv")
BULK_CIRRHOSIS <- file.path(INT_DIR, "dream_results_stage_cirrhosis.csv")

MIN_DONORS_PER_LR <- 30
MIN_DONORS_PER_STRATUM <- 10
TOP_N <- 100

# ----------------------------------------------------------------------------
# 0. Consolidate per-donor parquets -> merged TSV (avoids arrow R dep)
# ----------------------------------------------------------------------------
if (!file.exists(MERGED_TSV)) {
  cat(sprintf("[merge] consolidating parquets in %s\n", PER_DONOR_DIR))
  parquets <- list.files(PER_DONOR_DIR, pattern = "_lr_scores\\.parquet$",
                         full.names = TRUE)
  if (length(parquets) == 0) {
    stop("No per-donor parquets found in ", PER_DONOR_DIR)
  }
  cat(sprintf("[merge] %d parquets to consolidate\n", length(parquets)))
  # Run the consolidator via micromamba python (avoids R-side arrow dep)
  merge_py <- tempfile(fileext = ".py")
  writeLines(c(
    "import os, sys",
    "import pandas as pd",
    "from pathlib import Path",
    sprintf("PER_DONOR = Path(%s)", shQuote(PER_DONOR_DIR)),
    sprintf("OUT = Path(%s)",        shQuote(MERGED_TSV)),
    "parquets = sorted(PER_DONOR.glob('*_lr_scores.parquet'))",
    "dfs = [pd.read_parquet(p) for p in parquets]",
    "merged = pd.concat(dfs, ignore_index=True)",
    "print(f'shape={merged.shape}')",
    "merged.to_csv(OUT, sep='\\t', index=False, compression='gzip')"
  ), merge_py)
  cmd <- sprintf(
    "micromamba run -n rapids_singlecell python -u %s",
    shQuote(merge_py)
  )
  res <- system(cmd, intern = FALSE)
  if (res != 0 || !file.exists(MERGED_TSV)) {
    stop("Merge step failed (code=", res, ")")
  }
}

lr_long <- fread(MERGED_TSV)
cat(sprintf("[input] %d donor x LR rows from merged v2 TSV\n", nrow(lr_long)))

meta <- fread(META_V2)
cat(sprintf("[input] %d donors in v2 metadata\n", nrow(meta)))

covar_cols <- intersect(
  c("sample", "dataset", "disease_stage_coarse", "disease_stage_numeric",
    "F_stage_documented", "F_stage_inferred_v2", "F_stage_augmented_v2",
    "F_stage_source_v2", "macrophage_pseudotime_mean", "age", "sex_numeric"),
  names(meta)
)
lr_long <- merge(lr_long, meta[, ..covar_cols], by = "sample", all.x = TRUE)

lr_long[, score := -log10(pmax(magnitude_rank, 1e-4))]
stage_levels <- c("Healthy", "Steatosis", "Steatohepatitis", "Cirrhosis")
lr_long[, disease_stage_coarse := factor(disease_stage_coarse, levels = stage_levels)]
lr_long <- lr_long[!is.na(disease_stage_coarse)]
lr_long[, lr_pair := paste(ligand_complex, receptor_complex, sep = "__")]
lr_long[, ct_pair := paste(source, target, sep = "->")]

# ----------------------------------------------------------------------------
# Helper: fit one model
# ----------------------------------------------------------------------------
fit_one <- function(d, axis_term) {
  if (nrow(d) < MIN_DONORS_PER_LR) return(NULL)
  fixed <- c(axis_term)
  if (length(unique(d$dataset)) > 1) fixed <- c(fixed, "dataset")
  fixed <- c(fixed,
             "log10(pmax(n_source_cells,1))",
             "log10(pmax(n_target_cells,1))")
  if ("age" %in% names(d) &&
      length(unique(stats::na.omit(d$age))) > 1 &&
      sum(!is.na(d$age)) >= 10) {
    fixed <- c(fixed, "age")
  }
  if ("sex_numeric" %in% names(d) &&
      length(unique(stats::na.omit(d$sex_numeric))) > 1 &&
      sum(!is.na(d$sex_numeric)) >= 10) {
    fixed <- c(fixed, "sex_numeric")
  }
  rhs <- paste(fixed, collapse = " + ")
  use_random <- length(unique(d$dataset)) >= 3
  if (use_random) {
    rhs <- paste(rhs, "+ (1 | dataset)")
    fit <- try(lmerTest::lmer(stats::as.formula(paste("score ~", rhs)),
                              data = d, REML = TRUE), silent = TRUE)
  } else {
    fit <- try(stats::lm(stats::as.formula(paste("score ~", rhs)),
                         data = d), silent = TRUE)
  }
  if (inherits(fit, "try-error")) return(NULL)
  co <- try(summary(fit)$coefficients, silent = TRUE)
  if (inherits(co, "try-error")) return(NULL)
  rn <- rownames(co)
  hits <- rn[startsWith(rn, axis_term)]
  if (!length(hits)) return(NULL)
  out <- as.data.table(co[hits, , drop = FALSE], keep.rownames = "term")
  ncol_stats <- ncol(out) - 1L
  if (ncol_stats == 5L) {
    setnames(out, c("term", "Estimate", "StdErr", "df", "tval", "pval"))
  } else if (ncol_stats == 4L) {
    setnames(out, c("term", "Estimate", "StdErr", "tval", "pval"))
    out[, df := NA_real_]
  } else {
    return(NULL)
  }
  out[, n_donors := nrow(d)]
  out
}

fit_axis <- function(axis_term, only_with_value = NULL,
                     min_donors_per_stratum = MIN_DONORS_PER_STRATUM) {
  cat(sprintf("\n[axis] %s\n", axis_term))
  d_in <- lr_long
  if (!is.null(only_with_value)) {
    d_in <- d_in[!is.na(get(only_with_value))]
  }
  if (nrow(d_in) == 0) return(data.table())
  splits <- split(d_in, by = c("ct_pair", "lr_pair"), drop = TRUE)
  cat(sprintf("[axis] %d (ct_pair, lr_pair) groups\n", length(splits)))

  wrapper <- function(d) {
    if (axis_term == "disease_stage_coarse") {
      tab <- table(d$disease_stage_coarse)
      if (any(tab < min_donors_per_stratum)) return(NULL)
    }
    fit_res <- fit_one(d, axis_term)
    if (is.null(fit_res)) return(NULL)
    fit_res[, ct_pair := d$ct_pair[1]]
    fit_res[, lr_pair := d$lr_pair[1]]
    fit_res[, source := d$source[1]]
    fit_res[, target := d$target[1]]
    fit_res[, ligand_complex := d$ligand_complex[1]]
    fit_res[, receptor_complex := d$receptor_complex[1]]
    fit_res
  }

  t0 <- Sys.time()
  if (N_CORES > 1 && .Platform$OS.type != "windows") {
    res_list <- parallel::mclapply(splits, wrapper, mc.cores = N_CORES)
  } else {
    res_list <- lapply(splits, wrapper)
  }
  cat(sprintf("[axis] %s took %.1f min\n",
              axis_term, as.numeric(Sys.time() - t0, units = "mins")))
  out <- rbindlist(res_list[!sapply(res_list, is.null)], fill = TRUE)
  if (nrow(out) == 0) return(out)
  out[, padj_within_ct := stats::p.adjust(pval, method = "BH"),
      by = .(ct_pair, term)]
  out
}

# ============================================================================
# (a) 346 equivalent: stage-LR LMM (coarse / fstage / continuous)
# ============================================================================
cat("\n========== (346 v2): mixed-effects fits ==========\n")
res_coarse <- fit_axis("disease_stage_coarse")
fwrite(res_coarse, file.path(V2_STAGE_DIR, "stage_lr_lmm_coarse_v2.tsv"), sep = "\t")
cat(sprintf("[output] %d rows -> stage_lr_lmm_coarse_v2.tsv\n", nrow(res_coarse)))

# pick the best F-stage column for the v2 axis
fstage_col <- NULL
if ("F_stage_augmented_v2" %in% names(lr_long) &&
    sum(!is.na(lr_long$F_stage_augmented_v2)) >= 30) {
  fstage_col <- "F_stage_augmented_v2"
} else if ("F_stage_inferred_v2" %in% names(lr_long) &&
           sum(!is.na(lr_long$F_stage_inferred_v2)) >= 30) {
  fstage_col <- "F_stage_inferred_v2"
} else if ("F_stage_documented" %in% names(lr_long) &&
           sum(!is.na(lr_long$F_stage_documented)) >= 30) {
  fstage_col <- "F_stage_documented"
}
if (!is.null(fstage_col)) {
  cat(sprintf("[axis] F-stage column = %s (%d non-NA donors)\n",
              fstage_col, sum(!is.na(lr_long[[fstage_col]]))))
  lr_long[, F_stage_numeric := as.numeric(get(fstage_col))]
  res_fstage <- fit_axis("F_stage_numeric", only_with_value = "F_stage_numeric")
  fwrite(res_fstage, file.path(V2_STAGE_DIR, "stage_lr_lmm_fstage_v2.tsv"), sep = "\t")
  cat(sprintf("[output] %d rows -> stage_lr_lmm_fstage_v2.tsv\n", nrow(res_fstage)))
} else {
  cat("[axis] no F-stage column with sufficient coverage; skipping\n")
  res_fstage <- data.table()
}

if ("macrophage_pseudotime_mean" %in% names(lr_long) &&
    sum(!is.na(lr_long$macrophage_pseudotime_mean)) >= 30) {
  res_cont <- fit_axis("macrophage_pseudotime_mean",
                       only_with_value = "macrophage_pseudotime_mean")
} else {
  res_cont <- data.table()
}
fwrite(res_cont, file.path(V2_STAGE_DIR, "stage_lr_lmm_continuous_v2.tsv"), sep = "\t")
cat(sprintf("[output] %d rows -> stage_lr_lmm_continuous_v2.tsv\n", nrow(res_cont)))

# ============================================================================
# (b) 347 equivalent: leave-one-cohort-out replication
# ============================================================================
cat("\n========== (347 v2): LOO dataset replication ==========\n")
sh_term <- "disease_stage_coarseSteatohepatitis"
top <- res_coarse[term == sh_term][order(pval)][, head(.SD, TOP_N)]
cat(sprintf("[loo] selected top %d SH-vs-Healthy hits\n", nrow(top)))

datasets <- sort(unique(lr_long$dataset))
cat(sprintf("[loo] %d datasets: %s\n", length(datasets), paste(datasets, collapse = ", ")))

fit_one_holdout <- function(d, ds_held_out) {
  d <- d[dataset != ds_held_out]
  if (nrow(d) < 20) return(NULL)
  if (length(unique(d$disease_stage_coarse)) < 2) return(NULL)
  use_random <- length(unique(d$dataset)) >= 3
  fixed <- c("disease_stage_coarse",
             "log10(pmax(n_source_cells,1))",
             "log10(pmax(n_target_cells,1))")
  if (length(unique(d$dataset)) > 1) fixed <- c(fixed, "dataset")
  rhs <- paste(fixed, collapse = " + ")
  if (use_random) {
    fit <- try(lmerTest::lmer(stats::as.formula(paste("score ~", rhs, "+ (1|dataset)")),
                              data = d, REML = TRUE), silent = TRUE)
  } else {
    fit <- try(stats::lm(stats::as.formula(paste("score ~", rhs)), data = d),
               silent = TRUE)
  }
  if (inherits(fit, "try-error")) return(NULL)
  co <- try(summary(fit)$coefficients, silent = TRUE)
  if (inherits(co, "try-error")) return(NULL)
  rn <- rownames(co)
  hit_row <- rn[startsWith(rn, sh_term)]
  if (!length(hit_row)) return(NULL)
  list(estimate = co[hit_row, 1], pval = co[hit_row, ncol(co)])
}

run_loo_rows <- function() {
  out_rows <- list()
  for (ix in seq_len(nrow(top))) {
    hit <- top[ix]
    d_sub <- lr_long[ct_pair == hit$ct_pair & lr_pair == hit$lr_pair]
    if (nrow(d_sub) == 0) next
    for (ds in datasets) {
      f <- fit_one_holdout(d_sub, ds)
      if (is.null(f)) {
        replicates <- NA
        estimate <- NA_real_
        pval <- NA_real_
      } else {
        replicates <- (sign(f$estimate) == sign(hit$Estimate)) && (f$pval < 0.10)
        estimate <- f$estimate
        pval <- f$pval
      }
      out_rows[[length(out_rows) + 1]] <- data.table(
        ct_pair = hit$ct_pair, lr_pair = hit$lr_pair,
        full_estimate = hit$Estimate,
        full_padj_within_ct = hit$padj_within_ct,
        held_out = ds, holdout_estimate = estimate,
        holdout_pval = pval, replication_concordant = replicates
      )
    }
    if (ix %% 10 == 0) cat(sprintf("[loo] %d / %d\n", ix, nrow(top)))
  }
  rbindlist(out_rows)
}

if (nrow(top) > 0 && length(datasets) > 1) {
  loo_dt <- run_loo_rows()
  rep_rate <- loo_dt[, .(
    n_holdouts_total = length(datasets),
    n_holdouts_run = sum(!is.na(replication_concordant)),
    n_replicating = sum(replication_concordant, na.rm = TRUE),
    n_failed_fits = sum(is.na(replication_concordant))
  ), by = .(ct_pair, lr_pair)]
  rep_rate[, replication_rate_strict :=
                 n_replicating / n_holdouts_total]
  rep_rate[, replication_rate_among_run := ifelse(
    n_holdouts_run > 0, n_replicating / n_holdouts_run, NA_real_
  )]
  rep_rate[, replication_rate := replication_rate_among_run]
  fwrite(rep_rate,
         file.path(V2_STAGE_DIR, "loo_replication_rate_per_lr_v2.tsv"),
         sep = "\t")
  cat(sprintf("[output] %d LR pairs summarized -> loo_replication_rate_per_lr_v2.tsv\n",
              nrow(rep_rate)))
} else {
  rep_rate <- data.table()
  cat("[loo] skipped (insufficient top hits or single cohort)\n")
}

# ============================================================================
# (c) 348 equivalent: bulk concordance
# ============================================================================
cat("\n========== (348 v2): bulk concordance ==========\n")

load_bulk <- function(path) {
  if (!file.exists(path)) return(NULL)
  bulk <- fread(path)
  gene_col <- intersect(
    c("symbol", "gene_symbol", "hgnc_symbol", "gene", "Gene", "ID"),
    names(bulk)
  )[1]
  lfc_col <- intersect(c("logFC", "log2FoldChange", "lfc"), names(bulk))[1]
  padj_col <- intersect(c("padj", "adj.P.Val", "BH", "FDR"), names(bulk))[1]
  if (is.null(gene_col) || is.null(lfc_col) || is.null(padj_col)) return(NULL)
  out <- bulk[, .(gene = get(gene_col),
                  bulk_lfc = get(lfc_col),
                  bulk_padj = get(padj_col))]
  out <- out[!is.na(gene) & gene != ""]
  out <- out[order(bulk_padj)][, .SD[1], by = gene]
  out
}

bulk_all <- load_bulk(BULK_ALL)
bulk_st <- load_bulk(BULK_STEATOSIS)
bulk_sh <- load_bulk(BULK_SH)
bulk_cirr <- load_bulk(BULK_CIRRHOSIS)

score_concord <- function(lmm, anchor) {
  if (nrow(lmm) == 0 || is.null(anchor)) return(data.table())
  lmm <- copy(lmm)
  lmm[, ligand := sub("_.*", "", ligand_complex)]
  lmm[, receptor := sub("_.*", "", receptor_complex)]
  m_lig <- merge(lmm, anchor, by.x = "ligand", by.y = "gene", all.x = TRUE)
  setnames(m_lig, c("bulk_lfc", "bulk_padj"),
                  c("lig_bulk_lfc", "lig_bulk_padj"))
  m <- merge(m_lig, anchor, by.x = "receptor", by.y = "gene", all.x = TRUE)
  setnames(m, c("bulk_lfc", "bulk_padj"),
              c("rec_bulk_lfc", "rec_bulk_padj"))
  m[, sc_sign := sign(Estimate)]
  m[, lig_concordant := (sign(lig_bulk_lfc) == sc_sign) & (lig_bulk_padj < 0.1)]
  m[, rec_concordant := (sign(rec_bulk_lfc) == sc_sign) & (rec_bulk_padj < 0.1)]
  m[, both_concordant := lig_concordant & rec_concordant]
  m
}

score_save_axis <- function(term_label, anchor, axis_name) {
  slice <- res_coarse[term == term_label]
  if (nrow(slice) == 0) return(data.table())
  if (is.null(anchor)) {
    anchor <- bulk_all
    axis_tag <- paste0(axis_name, "_fallback_all_MASLD")
  } else {
    axis_tag <- axis_name
  }
  conc <- score_concord(slice, anchor)
  if (nrow(conc) > 0) conc[, axis := axis_tag]
  conc
}

steat_conc <- score_save_axis(
  "disease_stage_coarseSteatosis", bulk_st, "Steatosis_vs_Healthy"
)
sh_conc <- score_save_axis(
  "disease_stage_coarseSteatohepatitis", bulk_sh, "Steatohepatitis_vs_Healthy"
)
cirr_conc <- score_save_axis(
  "disease_stage_coarseCirrhosis", bulk_cirr, "Cirrhosis_vs_Healthy"
)

cont_conc <- data.table()
if (nrow(res_cont) > 0 && !is.null(bulk_all)) {
  cont_conc <- score_concord(res_cont, bulk_all)
  if (nrow(cont_conc) > 0) {
    cont_conc[, axis := "macrophage_pseudotime_all_MASLD"]
  }
}

fstage_conc <- data.table()
if (nrow(res_fstage) > 0 && !is.null(bulk_all)) {
  fstage_slice <- res_fstage[term == "F_stage_numeric"]
  fstage_conc <- score_concord(fstage_slice, bulk_all)
  if (nrow(fstage_conc) > 0) {
    fstage_conc[, axis := "F_stage_continuous_all_MASLD"]
  }
}

all_conc <- rbindlist(
  list(steat_conc, sh_conc, cirr_conc, cont_conc, fstage_conc),
  fill = TRUE
)
fwrite(all_conc, file.path(V2_STAGE_DIR, "lr_bulk_concordance_v2.tsv"), sep = "\t")
cat(sprintf("[output] %d total rows -> lr_bulk_concordance_v2.tsv\n", nrow(all_conc)))

if (nrow(all_conc) > 0) {
  hdr <- all_conc[!is.na(both_concordant),
                  .(n_testable = .N,
                    n_both = sum(both_concordant, na.rm = TRUE),
                    pct_both = 100 * mean(both_concordant, na.rm = TRUE),
                    pct_lig = 100 * mean(lig_concordant, na.rm = TRUE),
                    pct_rec = 100 * mean(rec_concordant, na.rm = TRUE)),
                  by = axis]
  cat("\n[headline] both_concordant per axis:\n")
  print(hdr)
}

# ============================================================================
# (d) 349 equivalent: figure (heatmap of top stage-progressive LR pairs)
# ============================================================================
cat("\n========== (349 v2): figure ==========\n")
if (nrow(res_coarse) > 0) {
  fig_top <- res_coarse[term == sh_term]
  if (nrow(fig_top) > 0) {
    fig_top[, neglog10p := -log10(pmax(pval, 1e-30))]
    fig_top[, rank_score := abs(Estimate) * neglog10p]
    fig_top <- fig_top[order(-rank_score)][, head(.SD, 60)]

    lr_long[, lr_pair := paste(ligand_complex, receptor_complex, sep = "__")]
    lr_long[, ct_pair := paste(source, target, sep = "->")]
    top_keys <- unique(fig_top[, .(ct_pair, lr_pair)])
    hm <- merge(lr_long, top_keys, by = c("ct_pair", "lr_pair"))
    hm_mean <- hm[, .(score_mean = mean(score, na.rm = TRUE),
                      n_donors = .N),
                  by = .(ct_pair, lr_pair, disease_stage_coarse)]
    hm_mean[, z := scale(score_mean)[, 1], by = .(ct_pair, lr_pair)]
    hm_mean[, disease_stage_coarse := factor(
      disease_stage_coarse, levels = stage_levels)]
    fig_top[, pair_id := paste(ct_pair, lr_pair, sep = " | ")]
    hm_mean[, pair_id := paste(ct_pair, lr_pair, sep = " | ")]
    ord <- fig_top[order(-Estimate)]$pair_id
    hm_mean[, pair_id := factor(pair_id, levels = ord)]

    p_heat <- ggplot(hm_mean, aes(x = disease_stage_coarse, y = pair_id, fill = z)) +
      geom_tile(color = "white", linewidth = 0.1) +
      scale_fill_gradient2(low = "#1565C0", mid = "#FFFFFF", high = "#C2185B",
                          midpoint = 0, name = "z (within LR)") +
      labs(x = NULL, y = NULL,
           title = "Top 60 stage-progressive LR pairs (v2 atlas)",
           subtitle = sprintf("Mean -log10(magnitude_rank) by stage (n_LR_total=%d)",
                              nrow(unique(hm[, .(ct_pair, lr_pair)])))) +
      theme_minimal(base_size = 7) +
      theme(axis.text.y = element_text(size = 5),
            axis.text.x = element_text(size = 7, angle = 30, hjust = 1),
            plot.title = element_text(face = "bold", size = 9),
            plot.subtitle = element_text(size = 7),
            legend.key.size = unit(0.3, "cm"))

    out_pdf <- file.path(FIG_DIR, "figS_stage_ccc_trajectory_v2_atlas.pdf")
    ggsave(out_pdf, p_heat, width = 7, height = 9, useDingbats = FALSE)
    cat(sprintf("[output] %s\n", out_pdf))
  } else {
    cat("[fig] no fig_top rows; skipping figure\n")
  }
} else {
  cat("[fig] no coarse fits; skipping figure\n")
}

cat("\n[done] 07_chain_346_349_v2 finished\n")
