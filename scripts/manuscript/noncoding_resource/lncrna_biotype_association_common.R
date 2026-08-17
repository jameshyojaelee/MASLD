#!/usr/bin/env Rscript

fit_biotype_association <- function(gene_table) {
  required <- c(
    "gene_id_versioned", "gene_type", "treat_positive", "AveExpr",
    "expression_variability", "gene_length_bp", "cohort_detection_count"
  )
  if (!all(required %in% names(gene_table))) {
    stop("Biotype-association fields are incomplete", call. = FALSE)
  }
  model_data <- gene_table[
    gene_type %in% c("protein_coding", "lncRNA") &
      is.finite(AveExpr) & is.finite(expression_variability) &
      is.finite(gene_length_bp) & is.finite(cohort_detection_count)
  ]
  if (nrow(model_data) < 100L || length(unique(model_data$treat_positive)) != 2L) {
    stop("Biotype-association model is not estimable", call. = FALSE)
  }
  model_data[, `:=`(
    is_lncrna = as.integer(gene_type == "lncRNA"),
    scaled_abundance = as.numeric(scale(AveExpr)),
    scaled_variability = as.numeric(scale(log1p(expression_variability))),
    scaled_log_length = as.numeric(scale(log(gene_length_bp))),
    scaled_cohort_detection = as.numeric(scale(cohort_detection_count))
  )]
  fit <- glm(
    treat_positive ~ is_lncrna + scaled_abundance + scaled_variability +
      scaled_log_length + scaled_cohort_detection,
    data = model_data,
    family = binomial()
  )
  coefficient <- summary(fit)$coefficients["is_lncrna", ]
  if (!all(is.finite(coefficient))) {
    stop("lncRNA coefficient is non-finite", call. = FALSE)
  }
  estimate <- unname(coefficient[["Estimate"]])
  standard_error <- unname(coefficient[["Std. Error"]])
  list(
    model_data = model_data,
    result = data.table::data.table(
      comparison = "lncRNA_vs_protein_coding",
      outcome = "all_gene_padj_below_0.05_and_abs_log2FC_above_0.50",
      n_genes = nrow(model_data),
      n_lncrna = sum(model_data$is_lncrna == 1L),
      n_protein_coding = sum(model_data$is_lncrna == 0L),
      log_odds_ratio = estimate,
      standard_error = standard_error,
      odds_ratio = exp(estimate),
      ci_lower = exp(estimate - 1.96 * standard_error),
      ci_upper = exp(estimate + 1.96 * standard_error),
      nominal_p = unname(coefficient[["Pr(>|z|)"]]),
      multiplicity_family = "one_prespecified_secondary_biotype_coefficient",
      interpretation = "descriptive_gene_level_association_not_biological_replication"
    )
  )
}
