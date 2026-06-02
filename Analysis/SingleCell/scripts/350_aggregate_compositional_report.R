#!/usr/bin/env Rscript
# S4 aggregate: build REPORT.md combining Fisher / scCODA / propeller p-values.
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

PROJECT <- Sys.getenv("MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUT     <- file.path(PROJECT, "Analysis/SingleCell/results_gpu_v2/proportion_compositional")
HEP     <- file.path(PROJECT, "Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes")

fisher_orig <- fread(file.path(HEP, "subtype_disease_enrichment.csv"))
fisher_orig[, subtype := as.character(subtype)]

read_opt <- function(p) {
  if (!file.exists(p)) { message("  missing: ", p); return(NULL) }
  fread(p)
}

sccoda <- read_opt(file.path(OUT, "sccoda_results.tsv"))
prop   <- read_opt(file.path(OUT, "propeller_results.tsv"))

merged_rows <- list()

# meta-subtype table
if (!is.null(sccoda)) {
  sc_meta <- sccoda[claim == "meta_subtype_Healthy_vs_MASLD"]
  if (!is.null(prop)) {
    pr_meta <- prop[claim == "meta_subtype_Healthy_vs_MASLD"]
  } else pr_meta <- data.table()
  if (nrow(sc_meta)) {
    merged_meta <- merge(
      sc_meta[, .(label, sccoda_final_param, sccoda_inclusion_prob, sccoda_credible)],
      if (nrow(pr_meta)) pr_meta[, .(label = cell_type, propeller_p = P.Value, propeller_padj = FDR,
                                     donor_asin_variance)] else data.table(label = character()),
      by = "label", all = TRUE
    )
    merged_meta[, level := "meta_subtype"]
    merged_rows[[length(merged_rows) + 1]] <- merged_meta
  }
}

if (!is.null(sccoda)) {
  sc_raw <- sccoda[grepl("^raw_subtype_", claim)]
  if (!is.null(prop)) pr_raw <- prop[claim == "raw_subtype_Healthy_vs_MASLD"] else pr_raw <- data.table()
  if (nrow(sc_raw)) {
    merged_raw <- merge(
      sc_raw[, .(label = as.character(label), sccoda_final_param, sccoda_inclusion_prob, sccoda_credible)],
      if (nrow(pr_raw)) pr_raw[, .(label = as.character(cell_type), propeller_p = P.Value, propeller_padj = FDR,
                                   donor_asin_variance)] else data.table(label = character()),
      by = "label", all = TRUE
    )
    merged_raw[, p_original_fisher := fisher_orig$fisher_pval[match(label, fisher_orig$subtype)]]
    merged_raw[, p_original_padj   := fisher_orig$fisher_padj[match(label, fisher_orig$subtype)]]
    merged_raw[, level := "raw_subtype"]
    merged_rows[[length(merged_rows) + 1]] <- merged_raw
  }
}

if (length(merged_rows)) {
  big <- rbindlist(merged_rows, use.names = TRUE, fill = TRUE)
  fwrite(big, file.path(OUT, "summary_table.tsv"), sep = "\t")
} else {
  big <- data.table()
}

# ---- REPORT.md ----
lines <- c(
  "# S4 — Compositional cell-type proportion retest",
  "",
  sprintf("Generated: %s", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  "",
  "## Methods",
  "",
  "Per-donor cell counts were modelled with two complementary methods:",
  "",
  "- **scCODA** (Büttner 2021): Bayesian Dirichlet-multinomial regression on per-donor count vectors; donor overdispersion is built into the likelihood. We use the binary Healthy-vs-MASLD covariate, 20,000 HMC iterations, and the most-abundant ubiquitous label as the reference.",
  "- **propeller** (Phipson 2022): empirical-Bayes moderated linear model on arcsin-sqrt-transformed per-donor proportions; robust = TRUE; donor is the replicate. Donor-level variance is reported as the arcsin-sqrt variance across donors.",
  "",
  "Original p-values are per-cluster Fisher exact tests on pooled cell counts (no donor structure) from `subtype_disease_enrichment.csv`.",
  "",
  "## Verdict",
  ""
)

if (nrow(big) > 0) {
  # Meta-subtype headline
  meta_tab <- big[level == "meta_subtype"]
  if (nrow(meta_tab)) {
    n_credible <- sum(meta_tab$sccoda_credible == TRUE, na.rm = TRUE)
    n_propeller_sig <- sum(meta_tab$propeller_padj < 0.05, na.rm = TRUE)
    lines <- c(lines,
      "### Meta-subtype (Progressor / Moderate / Stable / Neutral / Healthy / Disease-Neutral)",
      "",
      sprintf("- scCODA credible (inclusion-prob >= 0.95): **%d / %d** meta-labels", n_credible, nrow(meta_tab)),
      sprintf("- propeller FDR < 0.05: **%d / %d** meta-labels", n_propeller_sig, nrow(meta_tab)),
      "",
      "Per-label table:", ""
    )
    md_tab <- meta_tab
    md_tab[, sccoda_final_param := round(sccoda_final_param, 3)]
    md_tab[, sccoda_inclusion_prob := round(sccoda_inclusion_prob, 3)]
    if ("propeller_p" %in% names(md_tab)) {
      md_tab[, propeller_p := signif(propeller_p, 3)]
      md_tab[, propeller_padj := signif(propeller_padj, 3)]
    }
    lines <- c(lines, knitr::kable(md_tab, format = "pipe"), "")
  }

  raw_tab <- big[level == "raw_subtype"]
  if (nrow(raw_tab)) {
    n_credible <- sum(raw_tab$sccoda_credible == TRUE, na.rm = TRUE)
    n_propeller_sig <- sum(raw_tab$propeller_padj < 0.05, na.rm = TRUE)
    n_fisher_sig <- sum(raw_tab$p_original_padj < 0.05, na.rm = TRUE)
    lines <- c(lines,
      "### Raw 43-subtype",
      "",
      sprintf("- Original Fisher padj < 0.05:           **%d / %d**", n_fisher_sig, nrow(raw_tab)),
      sprintf("- scCODA credible (inclusion >= 0.95):   **%d / %d**", n_credible, nrow(raw_tab)),
      sprintf("- propeller FDR < 0.05:                  **%d / %d**", n_propeller_sig, nrow(raw_tab)),
      "",
      "Full table written to `summary_table.tsv`.",
      ""
    )
  }
} else {
  lines <- c(lines, "_No scCODA / propeller results loaded; both upstream jobs likely failed. See logs._", "")
}

lines <- c(lines,
  "## Files",
  "",
  "- `counts_meta_subtype_per_donor.tsv` — donor x meta-subtype cell counts",
  "- `counts_raw_subtype_per_donor.tsv`  — donor x raw-subtype cell counts",
  "- `sccoda_results.tsv`                — scCODA effect estimates + inclusion probabilities",
  "- `sccoda_meta_subtype.tsv` / `sccoda_raw_subtype.tsv` — raw scCODA effect tables",
  "- `propeller_results.tsv`             — propeller arcsin-LM p-values + donor variance",
  "- `summary_table.tsv`                 — joined comparison table",
  ""
)

writeLines(lines, file.path(OUT, "REPORT.md"))
message("[done] wrote ", file.path(OUT, "REPORT.md"))
