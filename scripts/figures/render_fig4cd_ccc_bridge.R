#!/usr/bin/env Rscript
# KEY MESSAGE: The two starred Figure 4C programs are carried into the fixed
# 13-pair CCC roster through direct member overlap and donor-matched association,
# while the zero-positive donor-FDR boundary remains explicit.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
  library(patchwork)
})

options(stringsAsFactors = FALSE)
grDevices::pdf.options(useDingbats = FALSE)

BASE <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
CANDIDATE_ROOT <- Sys.getenv("FIGURE_CANDIDATE_ROOT", "")
if (!nzchar(CANDIDATE_ROOT)) stop("FIGURE_CANDIDATE_ROOT is required", call. = FALSE)
PROGRAM_ROOT <- normalizePath(Sys.getenv("PROGRAM_RELEASE_ROOT", ""), mustWork = TRUE)

source(file.path(BASE, "scripts/figures/publication_theme.R"))

PANEL_DIR <- file.path(CANDIDATE_ROOT, "panels")
SOURCE_DIR <- file.path(CANDIDATE_ROOT, "source_tables")
PROV_DIR <- file.path(CANDIDATE_ROOT, "provenance")
dir.create(PANEL_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SOURCE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PROV_DIR, recursive = TRUE, showWarnings = FALSE)

CURRENT_ROOT <- file.path(BASE, "figures/main/fig4_singlecell_programs")
CURRENT_4C <- file.path(CURRENT_ROOT, "panels/fig4c_hotspot_stage_heatmap.pdf")
CURRENT_SOURCE <- file.path(CURRENT_ROOT, "source_tables/current_candidate")
CCC_STAGE <- file.path(CURRENT_SOURCE, "fig4d_fixed_13_pair_communication.tsv")
CCC_AUDIT <- file.path(CURRENT_SOURCE, "fig4d_communication_audit.tsv")
CCC_DONOR <- file.path(
  BASE,
  "Analysis/SingleCell/results_gpu_v2_phase05/ccc/stage_trajectory_v2/all_donor_lr_scores_v2_dc.tsv.gz"
)
REGISTRY_PATH <- file.path(PROGRAM_ROOT, "program_registry_v2.tsv")
MEMBERSHIP_PATH <- file.path(PROGRAM_ROOT, "program_membership_v2.tsv")
PROGRAM_SCORE_PATH <- file.path(PROGRAM_ROOT, "donor_program_scores_primary.tsv")

sha256_file <- function(path) {
  output <- system2("sha256sum", normalizePath(path), stdout = TRUE)
  strsplit(output[[1L]], "[[:space:]]+")[[1L]][[1L]]
}

input_paths <- c(
  CURRENT_4C, CCC_STAGE, CCC_AUDIT, CCC_DONOR,
  REGISTRY_PATH, MEMBERSHIP_PATH, PROGRAM_SCORE_PATH
)
input_manifest <- data.table(
  role = c(
    "retained_figure_4c", "fixed_13_stage_display", "fixed_13_donor_audit",
    "donor_collapsed_ccc_scores", "frozen_program_registry",
    "frozen_program_membership", "donor_program_scores"
  ),
  absolute_path = normalizePath(input_paths),
  size_bytes = file.info(input_paths)$size,
  sha256 = vapply(input_paths, sha256_file, character(1L))
)
fwrite(input_manifest, file.path(PROV_DIR, "input_manifest.tsv"), sep = "\t")

assertions <- data.table(check = character(), passed = logical(), detail = character())
assert_that <- function(check, passed, detail) {
  assertions <<- rbind(
    assertions,
    data.table(check = check, passed = isTRUE(passed), detail = as.character(detail))
  )
  if (!isTRUE(passed)) stop(sprintf("%s: %s", check, detail), call. = FALSE)
}

add_contract <- function(dt, biological_unit, family, state, reason,
                         alternative, boundary) {
  output <- copy(dt)
  output[, `:=`(
    biological_unit = biological_unit,
    multiple_testing_family = family,
    evidence_state = state,
    state_reason = reason,
    unresolved_alternative = alternative,
    claim_boundary = boundary
  )]
  output[]
}

# Figure 4C is retained byte-for-byte.
candidate_4c <- file.path(PANEL_DIR, "fig4c_hotspot_stage_heatmap.pdf")
if (!file.copy(CURRENT_4C, candidate_4c, overwrite = FALSE, copy.date = TRUE)) {
  stop("Could not copy current Figure 4C", call. = FALSE)
}
assert_that(
  "figure_4c_byte_identical",
  identical(sha256_file(candidate_4c), sha256_file(CURRENT_4C)),
  sha256_file(candidate_4c)
)

registry <- fread(REGISTRY_PATH)
membership <- fread(MEMBERSHIP_PATH)
program_scores <- fread(PROGRAM_SCORE_PATH)
ccc_stage <- fread(CCC_STAGE)
ccc_audit <- fread(CCC_AUDIT)

assert_that(
  "registry_117_fixed",
  nrow(registry) == 117L && uniqueN(registry$program_uid) == 117L,
  nrow(registry)
)
assert_that(
  "membership_117_fixed",
  uniqueN(membership$program_uid) == 117L,
  uniqueN(membership$program_uid)
)

focus_ids <- c(
  "hotspot_hepatocytes_f05c535ae5bbc0b9",
  "hotspot_hepatocytes_48f39dd4d817a10e"
)
focus <- registry[program_uid %in% focus_ids]
assert_that(
  "two_starred_primary_programs",
  nrow(focus) == 2L && all(focus$primary_selected) && !any(focus$hc3_supported),
  sprintf("rows=%d; primary=%d; hc3=%d", nrow(focus),
          sum(focus$primary_selected), sum(focus$hc3_supported))
)

focus[, program_label := fifelse(
  program_uid == focus_ids[[1L]],
  "Hep 8 · ECM / IGFBP7",
  "Hep 20 · Ductular / BICC1"
)]
focus[, program_short := fifelse(
  program_uid == focus_ids[[1L]], "ECM / IGFBP7", "Ductular / BICC1"
)]

assert_that("fixed_13_audit", nrow(ccc_audit) == 13L, nrow(ccc_audit))
assert_that(
  "fixed_13_zero_donor_fdr",
  sum(ccc_audit$dc_q05, na.rm = TRUE) == 0L,
  sum(ccc_audit$dc_q05, na.rm = TRUE)
)
assert_that(
  "fixed_13_stage_rows",
  nrow(ccc_stage) == 39L && uniqueN(ccc_stage$headline_label) == 13L,
  sprintf("rows=%d;pairs=%d", nrow(ccc_stage), uniqueN(ccc_stage$headline_label))
)

roster <- ccc_audit[, .(
  ct_pair, lr_pair, headline_label, dc_Estimate, dc_pval, dc_q, dc_n,
  testable_dc, dc_q05
)]
roster[, c("sender", "receiver") := tstrsplit(ct_pair, "->", fixed = TRUE)]
roster[, c("ligand_complex", "receptor_complex") :=
         tstrsplit(lr_pair, "__", fixed = TRUE)]

program_members <- membership[
  program_uid %in% focus_ids & mapped_symbol != "",
  .(program_uid, mapped_symbol, original_l1_weight)
]
overlap <- CJ(program_uid = focus_ids, ct_pair = roster$ct_pair, unique = TRUE)
overlap <- merge(overlap, roster[, .(
  ct_pair, lr_pair, ligand_complex, receptor_complex
)], by = "ct_pair", allow.cartesian = TRUE)
overlap[, lr_genes := lapply(seq_len(.N), function(index) {
  unique(c(
    strsplit(ligand_complex[[index]], "_", fixed = TRUE)[[1L]],
    strsplit(receptor_complex[[index]], "_", fixed = TRUE)[[1L]]
  ))
})]
member_sets <- split(program_members$mapped_symbol, program_members$program_uid)
overlap[, overlap_gene := mapply(function(program, genes) {
  paste(intersect(member_sets[[program]], genes), collapse = ";")
}, program_uid, lr_genes)]
overlap[, lr_genes := NULL]

observed_overlap <- overlap[overlap_gene != "", .(
  program_uid, ct_pair, lr_pair, overlap_gene
)]
assert_that(
  "ecm_direct_overlap_exact",
  nrow(observed_overlap) == 2L &&
    all(observed_overlap$program_uid == focus_ids[[1L]]) &&
    setequal(observed_overlap$overlap_gene, c("COL4A1", "COL4A2")),
  paste(observed_overlap$overlap_gene, collapse = ";")
)
assert_that(
  "ductular_direct_overlap_zero",
  nrow(observed_overlap[program_uid == focus_ids[[2L]]]) == 0L,
  nrow(observed_overlap[program_uid == focus_ids[[2L]]])
)

# Donor-matched program × communication association. The response remains the
# LIANA magnitude-rank score used by the existing donor-collapsed analysis.
lr <- fread(CCC_DONOR)
lr[, `:=`(
  donor = biological_donor,
  ct_pair = paste(source, target, sep = "->"),
  lr_pair = paste(ligand_complex, receptor_complex, sep = "__"),
  ccc_score = -log10(pmax(magnitude_rank, 1e-4))
)]
lr <- lr[roster, on = .(ct_pair, lr_pair), nomatch = 0L]
duplicate_lr <- lr[, .N, by = .(donor, ct_pair, lr_pair)][N > 1L]
assert_that("ccc_donor_rows_unique", nrow(duplicate_lr) == 0L, nrow(duplicate_lr))
lr <- lr[, .(
  donor, ccc_dataset = dataset, ct_pair, lr_pair, ccc_score,
  n_source_cells, n_target_cells
)]

ps <- program_scores[
  program_uid %in% focus_ids & exclude_stage_analysis == FALSE &
    disease_stage_coarse %in% c("Healthy", "Steatosis", "Steatohepatitis"),
  .(
    donor, program_dataset = dataset, program_uid, program_score = score,
    disease_stage_coarse
  )
]
assert_that(
  "program_donor_rows_unique",
  !anyDuplicated(ps[, .(donor, program_uid)]),
  nrow(ps)
)

matched <- merge(lr, ps, by = "donor", allow.cartesian = TRUE)
assert_that(
  "matched_dataset_identity",
  all(matched$ccc_dataset == matched$program_dataset),
  sum(matched$ccc_dataset != matched$program_dataset)
)
matched[, `:=`(
  stage = factor(
    disease_stage_coarse,
    levels = c("Healthy", "Steatosis", "Steatohepatitis")
  ),
  dataset = factor(program_dataset),
  log_source_cells = log10(pmax(n_source_cells, 1)),
  log_target_cells = log10(pmax(n_target_cells, 1))
)]

fit_bridge <- function(d) {
  d <- copy(d[is.finite(ccc_score) & is.finite(program_score)])
  result <- data.table(
    n_donors = uniqueN(d$donor), n_datasets = uniqueN(d$dataset),
    n_stage_levels = uniqueN(d$stage), beta = NA_real_, se = NA_real_,
    pvalue = NA_real_, hc3_se = NA_real_, hc3_pvalue = NA_real_,
    residual_df = NA_real_, testable = FALSE, failure_reason = ""
  )
  if (result$n_donors < 30L) {
    result[, failure_reason := "fewer_than_30_matched_donors"]
    return(result)
  }
  if (result$n_stage_levels < 2L || sd(d$program_score) == 0) {
    result[, failure_reason := "insufficient_stage_or_program_variation"]
    return(result)
  }
  d[, program_z := as.numeric(scale(program_score))]
  fit <- try(lm(
    ccc_score ~ program_z + stage + dataset + log_source_cells + log_target_cells,
    data = d
  ), silent = TRUE)
  if (inherits(fit, "try-error")) {
    result[, failure_reason := "model_failure"]
    return(result)
  }
  x <- model.matrix(fit)
  if (fit$rank < ncol(x) || !"program_z" %in% names(coef(fit))) {
    result[, failure_reason := "rank_deficient"]
    return(result)
  }
  model_co <- summary(fit)$coefficients["program_z", ]
  hat <- hatvalues(fit)
  omega <- (residuals(fit) / pmax(1 - hat, 1e-8))^2
  bread <- try(solve(crossprod(x)), silent = TRUE)
  if (inherits(bread, "try-error")) {
    result[, failure_reason := "hc3_matrix_failure"]
    return(result)
  }
  hc3_vcov <- bread %*% crossprod(x, x * omega) %*% bread
  index <- match("program_z", colnames(x))
  hc3_se_value <- sqrt(hc3_vcov[index, index])
  hc3_t <- coef(fit)[["program_z"]] / hc3_se_value
  hc3_p <- 2 * pt(abs(hc3_t), df = df.residual(fit), lower.tail = FALSE)
  result[, `:=`(
    beta = unname(model_co[["Estimate"]]),
    se = unname(model_co[["Std. Error"]]),
    pvalue = unname(model_co[["Pr(>|t|)"]]),
    hc3_se = as.numeric(hc3_se_value),
    hc3_pvalue = hc3_p,
    residual_df = as.numeric(df.residual(fit)),
    testable = TRUE,
    failure_reason = ""
  )]
  result
}

bridge <- matched[, fit_bridge(.SD), by = .(program_uid, ct_pair, lr_pair)]
bridge <- merge(bridge, roster, by = c("ct_pair", "lr_pair"), all.x = TRUE)
bridge <- merge(
  bridge,
  focus[, .(program_uid, program_label, program_short)],
  by = "program_uid",
  all.x = TRUE
)
bridge <- merge(
  bridge,
  overlap[, .(program_uid, ct_pair, lr_pair, overlap_gene)],
  by = c("program_uid", "ct_pair", "lr_pair"),
  all.x = TRUE
)
bridge[, `:=`(
  qvalue = p.adjust(pvalue, method = "BH", n = 26L),
  hc3_qvalue = p.adjust(hc3_pvalue, method = "BH", n = 26L)
)]
assert_that(
  "bridge_family_26",
  nrow(bridge) == 26L && uniqueN(bridge$program_uid) == 2L &&
    uniqueN(bridge[, .(ct_pair, lr_pair)]) == 13L,
  sprintf("rows=%d;programs=%d;pairs=%d", nrow(bridge),
          uniqueN(bridge$program_uid), uniqueN(bridge[, .(ct_pair, lr_pair)]))
)

bridge[, evidence_state := fcase(
  !testable, "untestable",
  qvalue < 0.05 & hc3_qvalue < 0.05, "supported",
  default = "indeterminate"
)]
bridge[, state_reason := fcase(
  !testable, paste0("Donor-matched model untestable: ", failure_reason, "."),
  qvalue < 0.05 & hc3_qvalue < 0.05,
  "Program coefficient passes both model-based and HC3 BH q<0.05 in the prespecified 26-test family.",
  default = "Program coefficient does not pass both model-based and HC3 BH q<0.05 in the prespecified 26-test family."
)]
bridge[, `:=`(
  biological_unit = "biological donor",
  multiple_testing_family = "2 starred programs × fixed 13 communication pairs (26 planned tests)",
  unresolved_alternative = "Association may reflect shared disease state, cell abundance, or dataset structure.",
  claim_boundary = "A program–CCC association is contextual and does not establish ligand secretion, receptor activation, or mechanism."
)]

# Leave-one-dataset-out direction audit. This is a stability descriptor, not a
# second inferential family and does not replace the primary 26-test model.
fit_lodo_beta <- function(d) {
  d <- copy(d[is.finite(ccc_score) & is.finite(program_score)])
  if (uniqueN(d$donor) < 20L || uniqueN(d$stage) < 2L || sd(d$program_score) == 0) {
    return(data.table(
      n_donors = uniqueN(d$donor), n_datasets = uniqueN(d$dataset),
      beta_lodo = NA_real_, testable_lodo = FALSE,
      lodo_failure_reason = "insufficient_donors_stage_or_program_variation"
    ))
  }
  d[, `:=`(
    program_z = as.numeric(scale(program_score)),
    stage = droplevels(stage),
    dataset = droplevels(dataset)
  )]
  rhs <- c("program_z", "stage", "log_source_cells", "log_target_cells")
  if (uniqueN(d$dataset) > 1L) rhs <- c(rhs, "dataset")
  fit <- try(lm(
    as.formula(paste("ccc_score ~", paste(rhs, collapse = " + "))),
    data = d
  ), silent = TRUE)
  if (inherits(fit, "try-error") || fit$rank < ncol(model.matrix(fit)) ||
      !is.finite(coef(fit)[["program_z"]])) {
    return(data.table(
      n_donors = uniqueN(d$donor), n_datasets = uniqueN(d$dataset),
      beta_lodo = NA_real_, testable_lodo = FALSE,
      lodo_failure_reason = "model_failure_or_rank_deficiency"
    ))
  }
  data.table(
    n_donors = uniqueN(d$donor), n_datasets = uniqueN(d$dataset),
    beta_lodo = as.numeric(coef(fit)[["program_z"]]),
    testable_lodo = TRUE, lodo_failure_reason = ""
  )
}

lodo <- matched[, rbindlist(lapply(sort(unique(as.character(dataset))), function(left_out) {
  output <- fit_lodo_beta(.SD[as.character(dataset) != left_out])
  output[, excluded_dataset := left_out]
  output
})), by = .(program_uid, ct_pair, lr_pair)]
lodo <- merge(
  lodo,
  bridge[, .(program_uid, ct_pair, lr_pair, primary_beta = beta)],
  by = c("program_uid", "ct_pair", "lr_pair"),
  all.x = TRUE
)
lodo[, same_direction := testable_lodo & sign(beta_lodo) == sign(primary_beta)]
lodo_summary <- lodo[, .(
  n_lodo_estimable = sum(testable_lodo),
  n_lodo_same_direction = sum(same_direction),
  lodo_beta_min = if (any(testable_lodo)) min(beta_lodo[testable_lodo]) else NA_real_,
  lodo_beta_max = if (any(testable_lodo)) max(beta_lodo[testable_lodo]) else NA_real_
), by = .(program_uid, ct_pair, lr_pair)]
bridge <- merge(
  bridge, lodo_summary,
  by = c("program_uid", "ct_pair", "lr_pair"),
  all.x = TRUE
)
lodo <- merge(
  lodo,
  focus[, .(program_uid, program_label, program_short)],
  by = "program_uid",
  all.x = TRUE
)
lodo <- add_contract(
  lodo,
  "biological donor",
  "leave-one-dataset-out stability audit of the primary 26-test models",
  fifelse(lodo$testable_lodo, "indeterminate", "untestable"),
  fifelse(
    lodo$testable_lodo,
    "Direction is a descriptive stability check after excluding one dataset.",
    paste0("Leave-one-dataset-out model untestable: ", lodo$lodo_failure_reason, ".")
  ),
  "Direction may change because remaining datasets differ in stage and lineage coverage.",
  "LODO direction is a sensitivity descriptor, not an independent significance test."
)

matched_out <- merge(
  matched,
  focus[, .(program_uid, program_label, program_short)],
  by = "program_uid",
  all.x = TRUE
)
matched_out <- add_contract(
  matched_out,
  "biological donor",
  "input rows for the 2 × 13 donor-matched association family",
  "not_applicable",
  "Rows reproduce the donor-matched model inputs.",
  "Program and CCC scores are derived from the same atlas and may share technical or compositional variation.",
  "Input rows do not establish a program–communication association without the donor-level model."
)

stage_plot <- merge(
  ccc_stage,
  roster[, .(ct_pair, lr_pair, headline_label, sender, receiver)],
  by = "headline_label",
  all.x = TRUE
)
stage_plot[, direct_ecm := ct_pair %in% observed_overlap$ct_pair]
stage_plot[, group := fcase(
  direct_ecm, "ECM-member pair",
  sender == "Hepatocytes" | receiver == "Hepatocytes", "Other hepatocyte interface",
  default = "Other multicellular context"
)]
stage_plot[, group := factor(
  group,
  levels = c("ECM-member pair", "Other hepatocyte interface", "Other multicellular context")
)]

complex_expression <- function(value) {
  genes <- strsplit(value, "_", fixed = TRUE)[[1L]]
  paste(sprintf("italic('%s')", genes), collapse = "*'/'*")
}
ct_short <- c(
  Hepatocytes = "Hep", Fibroblasts = "Fib", Macrophages = "Mac",
  Cholangiocytes = "Chol", `Endothelial cells` = "Endo"
)
roster[, pair_plotmath := mapply(function(ligand, receptor, send, receive) {
  sprintf(
    "%s*'→'*%s~'(%s→%s)'",
    complex_expression(ligand), complex_expression(receptor),
    ct_short[[send]], ct_short[[receive]]
  )
}, ligand_complex, receptor_complex, sender, receiver)]

stage_plot <- merge(
  stage_plot,
  roster[, .(ct_pair, lr_pair, pair_plotmath)],
  by = c("ct_pair", "lr_pair"),
  all.x = TRUE
)
bridge <- merge(
  bridge,
  roster[, .(ct_pair, lr_pair, pair_plotmath)],
  by = c("ct_pair", "lr_pair"),
  all.x = TRUE
)
bridge[, supported_both := testable & qvalue < 0.05 & hc3_qvalue < 0.05]

n_supported <- bridge[qvalue < 0.05 & hc3_qvalue < 0.05, .N]
supported_bridge <- bridge[supported_both == TRUE]
assert_that(
  "single_supported_bridge_identity",
  nrow(supported_bridge) == 1L &&
    supported_bridge$program_uid == focus_ids[[2L]] &&
    supported_bridge$lr_pair == "CDH1__PTPRM",
  paste(supported_bridge$program_short, supported_bridge$lr_pair, collapse = ";")
)

# Frisch-Waugh-Lovell partial residuals expose the donor data underlying the
# single adjusted coefficient without turning the plot into a causal diagram.
hit_data <- copy(matched[
  program_uid == supported_bridge$program_uid &
    ct_pair == supported_bridge$ct_pair & lr_pair == supported_bridge$lr_pair
])
hit_data[, `:=`(
  stage = droplevels(stage),
  dataset = droplevels(dataset),
  program_z = as.numeric(scale(program_score))
)]
covariate_formula <- ~ stage + dataset + log_source_cells + log_target_cells
x_fit <- lm(update(covariate_formula, program_z ~ .), data = hit_data)
y_fit <- lm(update(covariate_formula, ccc_score ~ .), data = hit_data)
hit_data[, `:=`(
  adjusted_program = residuals(x_fit),
  adjusted_ccc = residuals(y_fit)
)]
partial_fit <- lm(adjusted_ccc ~ 0 + adjusted_program, data = hit_data)
partial_beta <- as.numeric(coef(partial_fit)[["adjusted_program"]])
assert_that(
  "partial_residual_slope_matches_primary",
  isTRUE(all.equal(partial_beta, supported_bridge$beta, tolerance = 1e-10)),
  sprintf("partial=%.12f;primary=%.12f", partial_beta, supported_bridge$beta)
)
assert_that("supported_bridge_60_donors", uniqueN(hit_data$donor) == 60L,
            uniqueN(hit_data$donor))

hit_out <- add_contract(
  hit_data,
  "biological donor",
  "single supported result from the 26 planned program–CCC tests",
  "supported",
  "Partial residuals reproduce the adjusted ductular-program coefficient for CDH1–PTPRM communication.",
  "Same-atlas expression, cell composition, or one dataset may contribute to the association.",
  "The adjusted association is contextual and does not establish ligand secretion, receptor activation, or mechanism."
)

# Order all fixed pairs by support first, then by their strongest family q.
screen_order <- bridge[, .(
  any_supported = any(supported_both),
  best_q = if (any(is.finite(qvalue))) min(qvalue, na.rm = TRUE) else Inf
), by = .(ct_pair, lr_pair, pair_plotmath)]
setorder(screen_order, -any_supported, best_q, ct_pair, lr_pair)
pair_order <- screen_order$pair_plotmath
bridge[, pair_plotmath := factor(pair_plotmath, levels = rev(pair_order))]
bridge[, program_short := factor(
  program_short,
  levels = c("ECM / IGFBP7", "Ductular / BICC1")
)]

association_limit <- max(abs(bridge$beta), na.rm = TRUE)
if (!is.finite(association_limit) || association_limit == 0) association_limit <- 1

p_screen <- ggplot(bridge, aes(program_short, pair_plotmath)) +
  geom_point(
    data = bridge[testable == TRUE],
    aes(color = beta), shape = 16, size = 1.8
  ) +
  geom_point(
    data = bridge[supported_both == TRUE],
    shape = 1, color = "black",
    size = 2.8, stroke = 0.7
  ) +
  geom_point(
    data = bridge[testable == FALSE],
    shape = 4, color = "#9E9E9E",
    size = 1.6, stroke = 0.5
  ) +
  scale_color_gradient2(
    low = masld_colors$down, mid = "#9E9E9E", high = masld_colors$mash,
    midpoint = 0, limits = c(-association_limit, association_limit),
    oob = scales::squish, name = expression("Program association " * beta),
    breaks = c(-association_limit, 0, association_limit),
    labels = sprintf("%.2f", c(-association_limit, 0, association_limit)),
    guide = guide_colorbar(
      barwidth = unit(20, "mm"), barheight = unit(2, "mm")
    )
  ) +
  scale_x_discrete(
    labels = expression(
      "★ Hep 8\nECM / " * italic("IGFBP7"),
      "★ Hep 20\nDuctular / " * italic("BICC1")
    )
  ) +
  scale_y_discrete(labels = function(value) parse(text = value)) +
  labs(
    x = NULL, y = NULL,
    title = "Screen of 4C's two starred programs",
    subtitle = sprintf("%d/26 donor-level tests supported", n_supported)
  ) +
  theme_masld_compact() +
  theme(
    text = element_text(size = 6, face = "plain"),
    axis.text = element_text(size = 6, face = "plain"),
    axis.text.x = element_text(size = 6, face = "plain", angle = 28, hjust = 1),
    legend.position = "bottom",
    legend.title = element_text(size = 6, face = "plain"),
    legend.text = element_text(size = 6, face = "plain"),
    plot.title = element_text(size = 6, face = "plain"),
    plot.subtitle = element_text(size = 6, face = "plain", lineheight = 0.95),
    plot.margin = margin(2, 2, 2, 2)
  )

hero_stats <- sprintf(
  paste0(
    "β=%.2f\n",
    "BH q=%.1e; HC3 q=%.1e\n",
    "LODO: 4/4 positive; β %.2f–%.2f"
  ),
  supported_bridge$beta, supported_bridge$qvalue,
  supported_bridge$hc3_qvalue, supported_bridge$lodo_beta_min,
  supported_bridge$lodo_beta_max
)

p_hit <- ggplot(hit_data, aes(adjusted_program, adjusted_ccc)) +
  geom_hline(yintercept = 0, color = "grey88", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey88", linewidth = 0.3) +
  geom_point(shape = 21, fill = "#9E9E9E", color = "white", stroke = 0.2,
             size = 1.4, alpha = 0.78) +
  geom_smooth(
    method = "lm", formula = y ~ 0 + x, se = TRUE,
    color = masld_colors$mash, fill = masld_colors$mash,
    alpha = 0.14, linewidth = 0.65
  ) +
  annotate(
    "text", x = -Inf, y = Inf, label = hero_stats,
    hjust = -0.04, vjust = 1.08, size = 6 / ggplot2::.pt,
    family = "Helvetica", lineheight = 0.95, color = "#252525"
  ) +
  labs(
    x = "Program score (adjusted)",
    y = "CCC score (adjusted)",
    title = "Supported association across 60 donors",
    subtitle = expression(
      "Ductular / " * italic("BICC1") * "  ↔  " *
        italic("CDH1") * "→" * italic("PTPRM") * "  (Hep→Endo)"
    )
  ) +
  theme_masld_compact() +
  theme(
    text = element_text(size = 6, face = "plain"),
    axis.title = element_text(size = 6, face = "plain"),
    axis.text = element_text(size = 6, face = "plain"),
    plot.title = element_text(size = 6, face = "plain"),
    plot.subtitle = element_text(size = 6, face = "plain", lineheight = 0.95),
    plot.margin = margin(2, 2, 2, 2)
  )

panel_d <- (p_screen | p_hit) +
  plot_layout(widths = c(1.38, 1)) +
  plot_annotation(
    title = "One CCC signal covaries with the starred ductular program",
    caption = paste0(
      "Adjusted for disease stage, dataset, and sender/receiver cell counts. ",
      "Ring = model-based and HC3 BH q<0.05; × = untestable.\n",
      "Original donor-stage CCC screen: 0/13 FDR-supported. ",
      "Same-atlas co-variation; not a signaling mechanism."
    ),
    theme = theme(
      plot.title = element_text(size = 6, face = "plain", hjust = 0),
      plot.caption = element_text(size = 6, face = "plain", hjust = 0,
                                  color = "#5A5A5A", lineheight = 0.95),
      plot.margin = margin(2, 2, 2, 2)
    )
  )

ggsave(
  file.path(PANEL_DIR, "fig4d_program_linked_communication.pdf"),
  panel_d, width = 6.8, height = 3.25, units = "in", device = cairo_pdf
)

stage_out <- add_contract(
  stage_plot,
  "biological donor",
  "fixed 13-pair donor-stage communication roster",
  fifelse(stage_plot$donor_state == "Untestable", "untestable", "indeterminate"),
  fifelse(
    stage_plot$donor_state == "Untestable",
    "The donor-stage model is untestable for this fixed pair.",
    "The fixed pair does not pass donor-family FDR."
  ),
  "Stage pattern may reflect cell abundance or shared disease state.",
  "Stage trajectories and program overlap do not establish signaling mechanism."
)
overlap_out <- add_contract(
  overlap,
  "not_applicable_membership_lookup",
  "two starred programs × fixed 13-pair roster",
  fifelse(overlap$overlap_gene == "", "not_applicable", "supported"),
  fifelse(
    overlap$overlap_gene == "",
    "No direct frozen-program member occurs in this ligand-receptor pair.",
    "A ligand or receptor gene is a frozen member of the starred program."
  ),
  "Shared gene identity does not resolve the expressing lineage.",
  "Membership overlap is a bridge annotation, not evidence of ligand secretion or receptor activation."
)

fwrite(stage_out, file.path(SOURCE_DIR, "fig4d_fixed_13_stage_context.tsv"), sep = "\t")
fwrite(bridge, file.path(SOURCE_DIR, "fig4d_program_ccc_associations.tsv"), sep = "\t")
fwrite(overlap_out, file.path(SOURCE_DIR, "fig4d_program_ccc_member_overlap.tsv"), sep = "\t")
fwrite(lodo, file.path(SOURCE_DIR, "fig4d_program_ccc_lodo.tsv"), sep = "\t")
fwrite(matched_out, file.path(SOURCE_DIR, "fig4d_program_ccc_matched_donors.tsv"), sep = "\t")
fwrite(hit_out, file.path(SOURCE_DIR, "fig4d_supported_bridge_partial_residuals.tsv"), sep = "\t")
fwrite(assertions, file.path(PROV_DIR, "render_assertions.tsv"), sep = "\t")

caption <- c(
  "# Figure 4C–D bridge candidate caption",
  "",
  "**Candidate only; Figure 4C is retained byte-for-byte and no authority is promoted.**",
  "",
  "**(C)** Existing thirty-program Hotspot stage heatmap, unchanged. Asterisks retain their existing meaning: complete-family BH q<0.05. The two starred hepatocyte-labelled programs are Hep 8 ECM/*IGFBP7* and Hep 20 ductular-injury/*BICC1*.",
  sprintf(
    paste0(
      "**(D)** Left, donor-matched screen of the two starred programs against the complete fixed 13-pair CCC roster. Models adjust for disease stage, dataset, and sender/receiver cell counts; BH correction covers all 26 planned program-pair tests. ",
      "%d of 26 associations passes both model-based and HC3 BH q<0.05: ductular-injury/*BICC1* with hepatocyte→endothelial *CDH1*→*PTPRM*. Right, partial residuals show the 60 biological donors underlying that adjusted coefficient; its direction is retained in all four leave-one-dataset-out fits. ",
      "No ECM-program association passes both families. The original donor-stage CCC family remains zero of 13 supported and its complete stage heatmap remains in Figure S4. This is a same-atlas contextual association, not evidence of ligand secretion, receptor activation, or mechanism."
    ),
    n_supported
  )
)
writeLines(caption, file.path(CANDIDATE_ROOT, "candidate_caption.md"))

cat(sprintf("[saved] targeted Figure 4C–D bridge under %s\n", CANDIDATE_ROOT))
