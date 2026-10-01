#!/usr/bin/env Rscript

# Candidate Figure 3 insert: lncRNAs are abundant disease-associated
# transcripts, but are neither preferentially enriched after adjustment nor
# generally required to preserve prespecified program structure.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
})

options(digits = 17, scipen = 999)
set.seed(20260813)
grDevices::pdf.options(useDingbats = FALSE)

project_root <- normalizePath(
  Sys.getenv("MASLD_PROJECT_ROOT", "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"),
  mustWork = TRUE
)
source(file.path(project_root, "scripts/figures/publication_theme.R"), local = TRUE)
source(file.path(project_root, "scripts/figures/load_figure_data.R"), local = TRUE)

fail <- function(...) stop(..., call. = FALSE)
write_tsv_atomic <- function(x, path) {
  if (file.exists(path)) fail("Refusing to overwrite: ", path)
  tmp <- paste0(path, ".tmp.", Sys.getpid())
  fwrite(x, tmp, sep = "\t", quote = FALSE, na = "NA")
  if (!file.rename(tmp, path)) fail("Could not promote: ", path)
}
sha256_file <- function(path) {
  value <- system2("sha256sum", normalizePath(path, mustWork = TRUE), stdout = TRUE)
  strsplit(value[[1L]], "[[:space:]]+")[[1L]][[1L]]
}
plain_theme <- function() {
  theme_masld(base_size = 6) + theme_pub() +
    theme(
      text = element_text(size = 6, face = "plain", color = "black"),
      plot.title = element_text(size = 6, face = "plain", color = "black"),
      plot.subtitle = element_text(size = 6, face = "plain", color = "black"),
      legend.text = element_text(size = 6, face = "plain", color = "black"),
      legend.title = element_text(size = 6, face = "plain", color = "black")
    )
}
smd <- function(treated, control) {
  denominator <- sqrt((stats::var(treated) + stats::var(control)) / 2)
  if (!is.finite(denominator) || denominator == 0) return(0)
  (mean(treated) - mean(control)) / denominator
}

nc_root <- file.path(
  project_root, "RNA-seq/results/noncoding_resource/candidates/noncoding-bulk-2026-08-11"
)
# FIG3_LNCRNA_BIOTYPE_DIR / FIG3_LNCRNA_BULK_DIR: outputs of the two noncoding
# producers run on another bulk fit (run_fig3_lncrna_corrected.sbatch). With
# FIG3_LNCRNA_EXPECT_ADOPTED=0 the adopted-release counts below are replaced by
# agreement between the counts, covariate and association tables.
biotype_dir <- Sys.getenv("FIG3_LNCRNA_BIOTYPE_DIR",
  file.path(nc_root, "biotype_association/run_canonical_gate_2026-08-12"))
bulk_dir <- Sys.getenv("FIG3_LNCRNA_BULK_DIR",
  file.path(nc_root, "bulk_disease/run_canonical_gate_2026-08-12"))
expect_adopted <- Sys.getenv("FIG3_LNCRNA_EXPECT_ADOPTED", "1") != "0"
counts_path <- file.path(biotype_dir, "bulk_deg_counts_by_biotype.tsv")
association_path <- file.path(biotype_dir, "biotype_deg_association.tsv")
covariates_path <- file.path(biotype_dir, "biotype_association_gene_covariates.tsv")
bulk_path <- file.path(bulk_dir, "bulk_lncrna_results.tsv")
input_paths <- c(
  counts = counts_path,
  association = association_path,
  covariates = covariates_path,
  bulk_lncrna = bulk_path
)
if (any(!file.exists(input_paths))) fail("Missing input: ", paste(input_paths[!file.exists(input_paths)], collapse = ", "))

counts <- fread(counts_path)
association <- fread(association_path)
covariates <- fread(covariates_path)
bulk <- fread(bulk_path)

n_high_confidence <- sum(bulk$high_confidence == TRUE, na.rm = TRUE)
census <- covariates[, .(n_tested = .N, n_treat_positive = sum(treat_positive),
                         n_treat_up = sum(treat_positive & direction == "up"),
                         n_treat_down = sum(treat_positive & direction == "down")),
                     by = display_biotype]
census <- census[match(counts$display_biotype, display_biotype)]
if (nrow(counts) != 3L || nrow(association) != 1L ||
    !isTRUE(all.equal(as.data.frame(census[, -"display_biotype"]),
                      as.data.frame(counts[, .(n_tested, n_treat_positive, n_treat_up, n_treat_down)]),
                      check.attributes = FALSE)) ||
    association$n_genes != sum(counts[display_biotype != "other", n_tested])) {
  fail("Biotype counts, covariates and association disagree")
}
if (expect_adopted) {
  if (sum(counts$n_tested) != 23370L ||
      sum(counts$n_treat_positive) != 1347L) fail("Canonical biotype census drift")
  if (!identical(
    as.integer(counts[display_biotype == "protein_coding", c(n_treat_up, n_treat_down)]),
    c(744L, 142L)
  )) fail("Protein-coding direction census drift")
  if (!identical(
    as.integer(counts[display_biotype == "lncRNA", c(n_treat_up, n_treat_down)]),
    c(220L, 182L)
  )) fail("lncRNA direction census drift")
  if (abs(association$odds_ratio - 0.438898759683367) > 1e-12) {
    fail("Adjusted biotype association drift")
  }
  if (n_high_confidence != 279L) fail("High-confidence lncRNA count drift")
}

# Deterministic 1:1 propensity matching among canonical positives. Matching is
# descriptive at the gene level and conditions on the same four covariates as
# the prespecified enrichment model.
positive <- covariates[
  treat_positive == TRUE & display_biotype %in% c("protein_coding", "lncRNA")
]
if (expect_adopted && (nrow(positive[display_biotype == "protein_coding"]) != 886L ||
    nrow(positive[display_biotype == "lncRNA"]) != 402L)) fail("Positive matching universe drift")
positive[, is_lncrna := as.integer(display_biotype == "lncRNA")]
positive[, log_variability := log1p(expression_variability)]
positive[, log_length := log1p(gene_length_bp)]
propensity_fit <- glm(
  is_lncrna ~ scale(AveExpr) + scale(log_variability) + scale(log_length) +
    factor(cohort_detection_count),
  data = positive,
  family = binomial()
)
positive[, propensity_logit := as.numeric(predict(propensity_fit, type = "link"))]
# FIG3_LNCRNA_MATCH_METHOD:
#   propensity (default, the 2026-08-12 method): 1:1 nearest propensity logit,
#     without replacement, caliper 0.2 SD of the logit.
#   caliper_search: the same matching, trying calipers of 0.20, 0.15, 0.10 and
#     0.05 SD in that order and keeping the first that meets the unchanged 0.10
#     balance limit on all four covariates. The choice uses covariates only; the
#     direction outcome is computed after it. Introduced 2026-09-25 because the
#     0.2-SD match failed the balance limit on the corrected-control fit.
match_method <- Sys.getenv("FIG3_LNCRNA_MATCH_METHOD", "propensity")
if (!match_method %in% c("propensity", "caliper_search")) {
  fail("Unknown FIG3_LNCRNA_MATCH_METHOD: ", match_method)
}
match_covariates <- c("AveExpr", "log_variability", "log_length", "cohort_detection_count")
lnc <- positive[is_lncrna == 1L]
coding <- positive[is_lncrna == 0L]

match_at <- function(caliper_sd) {
  caliper <- caliper_sd * stats::sd(positive$propensity_logit)
  lnc_ordered <- lnc[order(-propensity_logit, gene_id_versioned)]
  available_coding <- coding[order(propensity_logit, gene_id_versioned)]
  selected <- vector("list", nrow(lnc_ordered))
  n_selected <- 0L
  for (i in seq_len(nrow(lnc_ordered))) {
    distances <- abs(available_coding$propensity_logit - lnc_ordered$propensity_logit[[i]])
    nearest <- which.min(distances)
    if (length(nearest) && distances[[nearest]] <= caliper) {
      n_selected <- n_selected + 1L
      selected[[n_selected]] <- data.table(
        lncrna_id = lnc_ordered$gene_id_versioned[[i]],
        coding_id = available_coding$gene_id_versioned[[nearest]],
        lncrna_logit = lnc_ordered$propensity_logit[[i]],
        coding_logit = available_coding$propensity_logit[[nearest]],
        distance = distances[[nearest]]
      )
      available_coding <- available_coding[-nearest]
    }
  }
  pairs <- rbindlist(selected[seq_len(n_selected)])
  lnc_match <- lnc[match(pairs$lncrna_id, gene_id_versioned)]
  coding_match <- coding[match(pairs$coding_id, gene_id_versioned)]
  if (anyNA(lnc_match$gene_id_versioned) || anyNA(coding_match$gene_id_versioned)) fail("Matched ID join failed")
  balance <- rbindlist(lapply(match_covariates, function(variable) data.table(
    variable = variable,
    pre_match_smd = smd(lnc[[variable]], coding[[variable]]),
    post_match_smd = smd(lnc_match[[variable]], coding_match[[variable]])
  )))
  balance[, `:=`(caliper_sd = caliper_sd, n_pairs = nrow(pairs))]
  list(pairs = pairs, lnc_match = lnc_match, coding_match = coding_match, balance = balance)
}

caliper_sequence <- if (match_method == "caliper_search") c(0.20, 0.15, 0.10, 0.05) else 0.20
balance_trace <- list()
matched <- NULL
for (caliper_sd in caliper_sequence) {
  attempt <- match_at(caliper_sd)
  balance_trace[[length(balance_trace) + 1L]] <- attempt$balance
  if (max(abs(attempt$balance$post_match_smd)) <= 0.10) {
    matched <- attempt
    break
  }
}
balance_trace <- rbindlist(balance_trace)
if (is.null(matched)) {
  print(balance_trace)
  fail("Post-match imbalance exceeds 0.10 at every caliper tried: ",
       paste(caliper_sequence, collapse = ", "))
}
pairs <- matched$pairs
lnc_match <- matched$lnc_match
coding_match <- matched$coding_match
balance <- balance_trace
balance[, `:=`(match_method = match_method, selected = caliper_sd == matched$balance$caliper_sd[[1L]])]
if (expect_adopted && nrow(pairs) < 360L) fail("Too few matched pairs: ", nrow(pairs))
pairs[, `:=`(pair_id = seq_len(.N), match_method = match_method,
             caliper_sd = matched$balance$caliper_sd[[1L]])]

pairs[, `:=`(
  lncrna_direction = lnc_match$direction,
  coding_direction = coding_match$direction,
  lncrna_up = lnc_match$direction == "up",
  coding_up = coding_match$direction == "up",
  lncrna_AveExpr = lnc_match$AveExpr,
  coding_AveExpr = coding_match$AveExpr,
  lncrna_expression_variability = lnc_match$expression_variability,
  coding_expression_variability = coding_match$expression_variability,
  lncrna_gene_length_bp = lnc_match$gene_length_bp,
  coding_gene_length_bp = coding_match$gene_length_bp
)]
discordant_lncrna_up <- sum(pairs$lncrna_up & !pairs$coding_up)
discordant_coding_up <- sum(!pairs$lncrna_up & pairs$coding_up)
if (discordant_lncrna_up == 0L || discordant_coding_up == 0L) fail("Matched direction OR is undefined")
matched_or <- discordant_lncrna_up / discordant_coding_up
matched_se <- sqrt(1 / discordant_lncrna_up + 1 / discordant_coding_up)
matched_ci <- exp(log(matched_or) + c(-1, 1) * 1.96 * matched_se)
matched_p <- stats::binom.test(
  discordant_lncrna_up,
  discordant_lncrna_up + discordant_coding_up,
  p = 0.5,
  alternative = "two.sided"
)$p.value

# Panel 1: apparent unadjusted enrichment reverses after prespecified adjustment.
lnc_count <- counts[display_biotype == "lncRNA"]
coding_count <- counts[display_biotype == "protein_coding"]
unadjusted_test <- fisher.test(matrix(c(
  lnc_count$n_treat_positive,
  lnc_count$n_tested - lnc_count$n_treat_positive,
  coding_count$n_treat_positive,
  coding_count$n_tested - coding_count$n_treat_positive
), nrow = 2, byrow = TRUE))
association_plot <- data.table(
  model = factor(
    c("Unadjusted", "Adjusted for\nexpression and detection"),
    levels = c("Adjusted for\nexpression and detection", "Unadjusted")
  ),
  odds_ratio = c(unname(unadjusted_test$estimate), association$odds_ratio),
  ci_lower = c(unadjusted_test$conf.int[[1L]], association$ci_lower),
  ci_upper = c(unadjusted_test$conf.int[[2L]], association$ci_upper)
)
p_association <- ggplot(association_plot, aes(x = odds_ratio, y = model)) +
  geom_vline(xintercept = 1, color = "#9E9E9E", linewidth = 0.3) +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), orientation = "y", width = 0.12, linewidth = 0.35) +
  geom_point(aes(shape = model, fill = model), size = 1.8, color = "black", stroke = 0.35) +
  scale_x_log10(limits = c(0.30, 1.40), breaks = c(0.5, 1)) +
  scale_shape_manual(
    values = c("Adjusted for\nexpression and detection" = 21, "Unadjusted" = 22),
    guide = "none"
  ) +
  scale_fill_manual(
    values = c("Adjusted for\nexpression and detection" = masld_colors$down, "Unadjusted" = "white"),
    guide = "none"
  ) +
  labs(
    title = "lncRNA DEG enrichment",
    x = "Odds of being an RNA-seq DEG\n(lncRNA / protein-coding)",
    y = NULL
  ) +
  plain_theme() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = margin(3, 5, 3, 3))

# Panel 2: bounded direction analysis among the canonical-positive genes.
direction_plot <- rbindlist(list(
  counts[display_biotype == "protein_coding", .(
    biotype = "Protein-coding", direction = "Upregulated\nin MASLD", n = n_treat_up, total = n_treat_positive
  )],
  counts[display_biotype == "protein_coding", .(
    biotype = "Protein-coding", direction = "Downregulated\nin MASLD", n = n_treat_down, total = n_treat_positive
  )],
  counts[display_biotype == "lncRNA", .(
    biotype = "lncRNA", direction = "Upregulated\nin MASLD", n = n_treat_up, total = n_treat_positive
  )],
  counts[display_biotype == "lncRNA", .(
    biotype = "lncRNA", direction = "Downregulated\nin MASLD", n = n_treat_down, total = n_treat_positive
  )]
))
direction_plot[, percent := 100 * n / total]
direction_plot[, biotype := factor(biotype, levels = c("Protein-coding", "lncRNA"))]
direction_plot[, direction := factor(direction, levels = c("Downregulated\nin MASLD", "Upregulated\nin MASLD"))]
p_direction <- ggplot(direction_plot, aes(x = biotype, y = percent, fill = direction)) +
  geom_col(width = 0.62, color = "white", linewidth = 0.25) +
  scale_fill_manual(values = c(
    "Downregulated\nin MASLD" = masld_colors$down,
    "Upregulated\nin MASLD" = masld_colors$up
  )) +
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 50, 100), labels = function(x) paste0(x, "%")) +
  labs(
    title = "Direction of RNA-seq DEGs",
    x = NULL,
    y = "Share of RNA-seq DEGs",
    fill = NULL
  ) +
  plain_theme() +
  theme(
    legend.position = "top", legend.direction = "horizontal",
    legend.key.width = unit(0.28, "cm"), legend.key.height = unit(0.20, "cm"),
    plot.margin = margin(1, 3, 0, 3)
  )

matched_direction_plot <- data.table(
  estimate = matched_or,
  lower = matched_ci[[1L]],
  upper = matched_ci[[2L]],
  row = "Covariate-matched\ngene pairs"
)
p_matched_direction <- ggplot(matched_direction_plot, aes(x = estimate, y = row)) +
  geom_vline(xintercept = 1, color = "#9E9E9E", linewidth = 0.3) +
  geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = 0.16, linewidth = 0.35) +
  geom_point(shape = 21, fill = masld_colors$up, color = "black", size = 1.8, stroke = 0.35) +
  scale_x_log10(
    limits = c(0.10, 1.25), breaks = c(0.1, 0.3, 1),
    labels = c("0.1", "0.3", "1")
  ) +
  labs(x = "Odds of upregulation in MASLD\n(lncRNA / protein-coding)", y = NULL) +
  plain_theme() +
  theme(
    axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    plot.margin = margin(0, 3, 3, 3)
  )

summary_source <- rbindlist(list(
  data.table(
    section = "composition", metric = c("tested_lncrna_share_percent", "positive_lncrna_share_percent", "high_confidence_lncrna"),
    value = c(100 * lnc_count$n_tested / sum(counts$n_tested), 100 * lnc_count$n_treat_positive / sum(counts$n_treat_positive), n_high_confidence),
    lower = NA_real_, upper = NA_real_, n = c(sum(counts$n_tested), sum(counts$n_treat_positive), lnc_count$n_treat_positive),
    unit = c("percent", "percent", "genes")
  ),
  association_plot[, .(
    section = "positive_call_association", metric = as.character(model), value = odds_ratio,
    lower = ci_lower, upper = ci_upper, n = association$n_genes[[1L]], unit = "odds_ratio"
  )],
  direction_plot[, .(
    section = "canonical_positive_direction", metric = paste(as.character(biotype), as.character(direction), sep = "_"),
    value = percent, lower = NA_real_, upper = NA_real_, n = n, unit = "percent"
  )],
  data.table(
    section = "matched_direction", metric = c("conditional_odds_ratio", "exact_mcnemar_p", "discordant_lncrna_up", "discordant_coding_up"),
    value = c(matched_or, matched_p, discordant_lncrna_up, discordant_coding_up),
    lower = c(matched_ci[[1L]], NA, NA, NA), upper = c(matched_ci[[2L]], NA, NA, NA),
    n = nrow(pairs), unit = c("odds_ratio", "p_value", "pairs", "pairs")
  ),
  data.table()
), fill = TRUE)

out_root <- Sys.getenv("FIG3_LNCRNA_OUT", FIG3_BULK_DIR)
panel_dir <- file.path(out_root, "panels")
source_dir <- file.path(out_root, "source_tables", "fig3_lncrna_role")
manifest_dir <- file.path(out_root, "manifests")
provenance_dir <- file.path(out_root, "provenance")
for (path in c(panel_dir, source_dir, manifest_dir, provenance_dir)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}
outputs <- c(
  pdf = file.path(panel_dir, "fig3_lncrna_role_candidate.pdf"),
  summary = file.path(source_dir, "summary.tsv"),
  pairs = file.path(source_dir, "matched_pairs.tsv"),
  balance = file.path(source_dir, "match_balance.tsv"),
  session = file.path(provenance_dir, "fig3_lncrna_role_sessionInfo.txt"),
  manifest = file.path(manifest_dir, "fig3_lncrna_role_candidate_manifest.tsv")
)
if (any(file.exists(outputs))) fail("Refusing existing candidate output: ", paste(outputs[file.exists(outputs)], collapse = ", "))

pdf_tmp <- paste0(outputs[["pdf"]], ".tmp.", Sys.getpid())
grDevices::cairo_pdf(pdf_tmp, width = 4.25, height = 2.35, family = "Helvetica")
grid.newpage()
layout <- grid.layout(1, 2, widths = unit(c(1.05, 1.15), "null"))
pushViewport(viewport(layout = layout))
print(p_association, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
pushViewport(viewport(
  layout.pos.row = 1, layout.pos.col = 2,
  layout = grid.layout(2, 1, heights = unit(c(1.35, 0.65), "null"))
))
print(p_direction, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_matched_direction, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
upViewport()
dev.off()
if (!file.rename(pdf_tmp, outputs[["pdf"]])) fail("Could not promote PDF")

write_tsv_atomic(summary_source, outputs[["summary"]])
write_tsv_atomic(pairs, outputs[["pairs"]])
write_tsv_atomic(balance, outputs[["balance"]])
writeLines(capture.output(sessionInfo()), outputs[["session"]], useBytes = TRUE)

manifest <- data.table(
  role = c(names(input_paths), names(outputs)[names(outputs) != "manifest"]),
  direction = c(rep("input", length(input_paths)), rep("output", length(outputs) - 1L)),
  path = c(input_paths, outputs[names(outputs) != "manifest"])
)
manifest[, `:=`(size_bytes = file.info(path)$size, sha256 = vapply(path, sha256_file, character(1)))]
write_tsv_atomic(manifest, outputs[["manifest"]])

message(sprintf(
  "PASS: %d matched pairs (%s, caliper %.2f SD); matched direction OR %.3f (95%% CI %.3f-%.3f); max |SMD| %.4f",
  nrow(pairs), match_method, pairs$caliper_sd[[1L]], matched_or, matched_ci[[1L]], matched_ci[[2L]],
  max(abs(balance[selected == TRUE, post_match_smd]))
))
message("Wrote: ", outputs[["pdf"]])
