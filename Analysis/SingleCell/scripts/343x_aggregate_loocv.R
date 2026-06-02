#!/usr/bin/env Rscript
# 343x_aggregate_loocv.R  (S2 deliverable #3)
#
# Aggregate the 10 donor-jackknife + 7 cohort-out scVI retrains produced by
# 343v / 343w into an honest QWK distribution.
#
# Reads:
#   results_gpu_v2/scvi_validation/donor_jackknife/<donor>/summary.json
#   results_gpu_v2/scvi_validation/cohort_out/<cohort>/summary.json
#
# Writes:
#   results_gpu_v2/scvi_validation/donor_jackknife_summary.tsv
#   results_gpu_v2/scvi_validation/cohort_out_summary.tsv
#   results_gpu_v2/scvi_validation/REPORT.md
#
# Reports the honest QWK alongside the leaked QWK = 0.76 number from
# 343b kNN-on-fixed-latent so the comparison is unambiguous.

suppressPackageStartupMessages({
  library(jsonlite)
  library(dplyr)
})

WT_ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation"
OUT_DIR <- file.path(WT_ROOT, "Analysis/SingleCell/results_gpu_v2/scvi_validation")
DONOR_DIR  <- file.path(OUT_DIR, "donor_jackknife")
COHORT_DIR <- file.path(OUT_DIR, "cohort_out")

# ---- donor jackknife -------------------------------------------------------
donor_rows <- list()
if (dir.exists(DONOR_DIR)) {
  for (d in list.dirs(DONOR_DIR, recursive = FALSE)) {
    f <- file.path(d, "summary.json")
    if (file.exists(f)) {
      # Python json writes np.nan as literal `NaN` which is invalid JSON.
      # Sanitize before parsing.
      txt <- paste(readLines(f, warn = FALSE), collapse = "\n")
      txt <- gsub("\\bNaN\\b", "null", txt)
      txt <- gsub("\\bInfinity\\b", "null", txt)
      txt <- gsub("\\b-Infinity\\b", "null", txt)
      x <- jsonlite::fromJSON(txt)
      # Coerce list values to length-1 vectors / handle NULL → NA so as.data.frame
      # doesn't choke on heterogeneous summary.json files.
      x_norm <- lapply(x, function(v) {
        if (is.null(v) || length(v) == 0) NA
        else if (length(v) > 1) paste(unlist(v), collapse = ";")
        else v
      })
      donor_rows[[length(donor_rows) + 1]] <- tryCatch(
        as.data.frame(x_norm, stringsAsFactors = FALSE),
        error = function(e) { warning("skipping ", f, ": ", conditionMessage(e)); NULL }
      )
    }
  }
}
donor_df <- if (length(donor_rows)) dplyr::bind_rows(donor_rows) else data.frame()

if (nrow(donor_df) >= 2 && all(c("F_true", "F_pred") %in% colnames(donor_df))) {
  donor_df$F_true <- as.integer(donor_df$F_true)
  donor_df$F_pred <- as.integer(donor_df$F_pred)
  donor_df$abs_err <- abs(donor_df$F_true - donor_df$F_pred)
  # quadratic-weighted kappa (psych::cohen.kappa) — implement manually to avoid dep
  qwk_manual <- function(y, p, k = 5) {
    cm <- table(factor(y, levels = 0:(k-1)), factor(p, levels = 0:(k-1)))
    cm <- as.matrix(cm)
    n  <- sum(cm)
    w  <- outer(0:(k-1), 0:(k-1), function(i, j) ((i - j)^2) / ((k - 1)^2))
    Eh <- outer(rowSums(cm), colSums(cm)) / n
    1 - sum(w * cm) / sum(w * Eh)
  }
  donor_qwk <- tryCatch(qwk_manual(donor_df$F_true, donor_df$F_pred),
                        error = function(e) NA_real_)
} else {
  donor_qwk <- NA_real_
}

write.table(donor_df, file.path(OUT_DIR, "donor_jackknife_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- cohort-out ------------------------------------------------------------
cohort_rows <- list()
if (dir.exists(COHORT_DIR)) {
  for (d in list.dirs(COHORT_DIR, recursive = FALSE)) {
    f <- file.path(d, "summary.json")
    if (file.exists(f)) {
      # Python json writes np.nan as literal `NaN` which is invalid JSON.
      # Sanitize before parsing.
      txt <- paste(readLines(f, warn = FALSE), collapse = "\n")
      txt <- gsub("\\bNaN\\b", "null", txt)
      txt <- gsub("\\bInfinity\\b", "null", txt)
      txt <- gsub("\\b-Infinity\\b", "null", txt)
      x <- jsonlite::fromJSON(txt)
      x_norm <- lapply(x, function(v) {
        if (is.null(v) || length(v) == 0) NA
        else if (length(v) > 1) paste(unlist(v), collapse = ";")
        else v
      })
      cohort_rows[[length(cohort_rows) + 1]] <- tryCatch(
        as.data.frame(x_norm, stringsAsFactors = FALSE),
        error = function(e) { warning("skipping ", f, ": ", conditionMessage(e)); NULL }
      )
    }
  }
}
cohort_df <- if (length(cohort_rows)) dplyr::bind_rows(cohort_rows) else data.frame()
write.table(cohort_df, file.path(OUT_DIR, "cohort_out_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- REPORT.md -------------------------------------------------------------
md <- c(
  "# scVI proper-LOOCV REPORT (S2)",
  "",
  paste0("Date: ", format(Sys.Date())),
  "",
  "## Headline",
  "",
  sprintf("- Donor-jackknife honest QWK (n = %d retrains): **%s**",
          nrow(donor_df),
          if (is.na(donor_qwk)) "NA" else sprintf("%.3f", donor_qwk)),
  sprintf("- Cohort-out honest QWK (mean across %d cohorts with labels): **%s**",
          if (nrow(cohort_df)) sum(!is.na(cohort_df$qwk)) else 0L,
          if (nrow(cohort_df) && any(!is.na(cohort_df$qwk)))
            sprintf("%.3f", mean(cohort_df$qwk, na.rm = TRUE)) else "NA"),
  "- Leaked QWK from 343b (kNN-5 on fixed scVI latent, all donors in scVI training): 0.760",
  "",
  "## Honesty caveat",
  "",
  "Donor-jackknife and cohort-out each retrain scVI from scratch on the held-out-removed AnnData and run latent inference on the unseen cells via `prepare_query_anndata + get_latent_representation`. No parameter updates happen during projection.",
  "",
  "## Donor jackknife details",
  ""
)
if (nrow(donor_df)) {
  md <- c(md, "| sample | F_true | F_pred | abs_err | n_train_donors_labelled | n_held_cells |",
          "|---|---|---|---|---|---|")
  for (i in seq_len(nrow(donor_df))) {
    md <- c(md, sprintf("| %s | %s | %s | %s | %s | %s |",
                        donor_df$sample[i], donor_df$F_true[i], donor_df$F_pred[i],
                        donor_df$abs_err[i], donor_df$n_train_donors_labelled[i],
                        donor_df$n_held_cells[i]))
  }
} else {
  md <- c(md, "(no donor jackknife jobs completed yet)")
}

md <- c(md, "", "## Cohort-out details", "")
if (nrow(cohort_df)) {
  md <- c(md,
          "| cohort | QWK | latent_KL | marker_rho | n_held_donors_labelled | n_train_donors_labelled | n_held_cells |",
          "|---|---|---|---|---|---|---|")
  for (i in seq_len(nrow(cohort_df))) {
    md <- c(md, sprintf("| %s | %s | %.3f | %s | %s | %s | %s |",
                        cohort_df$hold_out_cohort[i],
                        if (is.na(cohort_df$qwk[i])) "NA" else sprintf("%.3f", cohort_df$qwk[i]),
                        cohort_df$latent_kl[i],
                        if (is.null(cohort_df$marker_preservation_rho[i]) || is.na(cohort_df$marker_preservation_rho[i])) "NA"
                          else sprintf("%.3f", as.numeric(cohort_df$marker_preservation_rho[i])),
                        cohort_df$n_held_donors_labelled[i],
                        cohort_df$n_train_donors_labelled[i],
                        cohort_df$n_held_cells[i]))
  }
} else {
  md <- c(md, "(no cohort-out jobs completed yet)")
}

md <- c(md, "",
        "## Interpretation guide",
        "",
        "- Honest QWK substantially below 0.76 ⇒ leaked latent inflated the prior estimate; report the honest value as the F-stage generalization benchmark and revise `memory/MEMORY.md`.",
        "- Per-cohort QWK with high variance ⇒ F-stage signal is cohort-specific (e.g. only generalizes within Andrews); flag this for the paper.",
        "- High latent_KL on a cohort with low QWK ⇒ scVI failed to integrate that cohort; the projection is unreliable, not the F-stage classifier.",
        "- Low marker_preservation_rho on a cohort ⇒ that cohort's hepatocyte transcriptome diverges from the training-cohort hepatocytes (could be technical, e.g. snRNA vs scRNA).")

writeLines(md, file.path(OUT_DIR, "REPORT.md"))
cat("wrote", file.path(OUT_DIR, "REPORT.md"), "\n")
