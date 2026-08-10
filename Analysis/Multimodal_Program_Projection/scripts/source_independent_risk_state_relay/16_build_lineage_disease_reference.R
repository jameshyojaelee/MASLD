#!/usr/bin/env Rscript
# Freeze the Plan 45 lineage-resolved human MASLD disease reference.
#
# This is a new analysis, not a reuse of the failed Plan 43/44 global axes or
# the older coarse-stage pseudobulk output. Raw run-keyed counts are collapsed
# to biological donors before modeling, and the >=50-cell gate uses the exact
# count for the lineage being analyzed.

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(jsonlite)
  library(limma)
})

project_root <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
candidate_id <- Sys.getenv(
  "PLAN45_DISEASE_REFERENCE_CANDIDATE_ID",
  unset = "source-independent-risk-state-relay-lineage-disease-reference-2026-08-10"
)
if (!grepl("^source-independent-risk-state-relay-[A-Za-z0-9._-]+$", candidate_id)) {
  stop("Invalid PLAN45_DISEASE_REFERENCE_CANDIDATE_ID")
}
candidate_parent <- file.path(
  project_root, "Analysis/Multimodal_Program_Projection/candidates"
)
candidate_root <- file.path(candidate_parent, candidate_id)
if (dir.exists(candidate_root)) stop("Refusing to overwrite candidate: ", candidate_root)

pb_root <- file.path(project_root, "Analysis/SingleCell/results_gpu_v2/pseudobulk")
metadata_path <- file.path(
  project_root,
  "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_metadata_extended.tsv"
)
collapse_path <- file.path(project_root, "Analysis/SingleCell/scripts/lib_donor_collapse.R")
annotation_path <- file.path(
  project_root, "Analysis/SingleCell/integration/scripts/annotate_celltypist.py"
)
lineage_reference_root <- file.path(
  candidate_parent,
  "source-independent-risk-state-relay-lineage-observability-2026-08-10"
)
lineage_reference_seal <- file.path(lineage_reference_root, "LINEAGE_REFERENCE_SEALED.json")
required_inputs <- c(metadata_path, collapse_path, annotation_path, lineage_reference_seal)
if (any(!file.exists(required_inputs))) {
  stop("Required lineage-disease input absent: ",
       paste(required_inputs[!file.exists(required_inputs)], collapse = "; "))
}
source(collapse_path)

lineage_spec <- data.table(
  cell_type = c(
    "Hepatocytes", "Macrophages", "Fibroblasts", "Cholangiocytes",
    "Endothelial_cells"
  ),
  cell_count_column = c(
    "n_Hepatocytes", "n_Macrophages", "n_Fibroblasts", "n_Cholangiocytes",
    "n_Endothelial_cells"
  ),
  recipient_role = c(
    "primary_organoid", "primary_organoid", "primary_organoid",
    "primary_organoid", "optional_recipient"
  )
)

# These genes were explicitly used to validate the original CellTypist lineage
# annotation. They are excluded from the primary state projection so a change
# in annotation identity cannot masquerade as a disease relay.
identity_markers <- unique(c(
  "ALB", "APOB", "APOC3", "CYP3A4",
  "KRT19", "EPCAM", "KRT7",
  "PECAM1", "VWF", "CDH5",
  "CD68", "MARCO", "CD163",
  "ACTA2", "COL1A1", "COL3A1"
))

primary_stages <- c("Healthy", "Steatosis", "Steatohepatitis")
source_uncertain_dataset <- "GSE189600"
primary_min_cells <- 50L
sensitivity_min_cells <- 20L
treat_lfc <- 0.25

sha256 <- function(path) {
  output <- system2("sha256sum", args = path, stdout = TRUE, stderr = TRUE)
  if (length(output) != 1L) stop("sha256sum failed for ", path)
  strsplit(output, "[[:space:]]+")[[1L]][[1L]]
}

first_nonmissing <- function(x) {
  y <- x[!is.na(x) & as.character(x) != ""]
  if (length(y)) y[[1L]] else NA
}

parse_flag <- function(x) {
  tolower(trimws(as.character(x))) %chin% c("true", "t", "1", "yes", "y")
}

as_numeric_zero <- function(x) {
  z <- suppressWarnings(as.numeric(x))
  z[!is.finite(z)] <- 0
  z
}

assert <- function(condition, ...) {
  if (!isTRUE(condition)) stop(..., call. = FALSE)
}

atomic_fwrite <- function(x, path) {
  tmp <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  fwrite(x, tmp, sep = "\t", quote = FALSE, na = "")
  if (!file.rename(tmp, path)) stop("Atomic rename failed for ", path)
}

bh_adjust <- function(p) {
  out <- rep(NA_real_, length(p))
  ok <- is.finite(p)
  out[ok] <- p.adjust(p[ok], method = "BH")
  out
}

build_design <- function(meta) {
  d <- copy(meta)
  if (uniqueN(d$dataset) > 1L) {
    d[, dataset := factor(dataset, levels = sort(unique(dataset)))]
    design <- model.matrix(~ dataset + disease_binary, data = d)
    residual_disease <- residuals(lm(disease_binary ~ dataset, data = d))
  } else {
    design <- model.matrix(~ disease_binary, data = d)
    residual_disease <- residuals(lm(disease_binary ~ 1, data = d))
  }
  coef_name <- "disease_binary"
  if (!coef_name %in% colnames(design)) stop("Disease coefficient absent from design")
  list(
    meta = d,
    design = design,
    coef_index = match(coef_name, colnames(design)),
    rank = qr(design)$rank,
    n_columns = ncol(design),
    condition_number = kappa(design, exact = TRUE),
    residual_disease_variance = if (length(residual_disease) > 1L) {
      var(residual_disease)
    } else {
      NA_real_
    }
  )
}

empty_effect_table <- function(genes, cell_type, model_id, reason) {
  data.table(
    cell_type = cell_type,
    model_id = model_id,
    gene_symbol = genes,
    filter_by_expr = FALSE,
    logFC = NA_real_, moderated_se = NA_real_, AveExpr = NA_real_,
    t = NA_real_, pvalue = NA_real_, qvalue = NA_real_, B = NA_real_,
    treat_lfc = treat_lfc, treat_t = NA_real_, treat_pvalue = NA_real_,
    treat_qvalue = NA_real_, estimable = FALSE, failure_reason = reason
  )
}

fit_voom <- function(counts, meta, cell_type, model_id,
                     min_total = 12L, min_each = 5L, min_residual_df = 5L) {
  genes <- rownames(counts)
  base_audit <- data.table(
    cell_type = cell_type,
    model_id = model_id,
    n_donors = nrow(meta),
    n_controls = sum(meta$disease_binary == 0L),
    n_disease = sum(meta$disease_binary == 1L),
    n_datasets = uniqueN(meta$dataset),
    n_mixed_datasets = meta[, uniqueN(disease_binary), by = dataset][V1 == 2L, .N],
    design_rank = NA_integer_, n_model_columns = NA_integer_,
    residual_df = NA_integer_, condition_number = NA_real_,
    residual_disease_variance = NA_real_, n_filter_by_expr = 0L,
    estimable = FALSE, failure_reason = ""
  )
  if (nrow(meta) < min_total || sum(meta$disease_binary == 0L) < min_each ||
      sum(meta$disease_binary == 1L) < min_each) {
    base_audit[, failure_reason := "insufficient_biological_donors"]
    return(list(
      effects = empty_effect_table(genes, cell_type, model_id,
                                   "insufficient_biological_donors"),
      audit = base_audit
    ))
  }
  if (base_audit$n_mixed_datasets < 1L) {
    base_audit[, failure_reason := "disease_not_identified_within_dataset"]
    return(list(
      effects = empty_effect_table(genes, cell_type, model_id,
                                   "disease_not_identified_within_dataset"),
      audit = base_audit
    ))
  }
  design_obj <- build_design(meta)
  design <- design_obj$design
  base_audit[, `:=`(
    design_rank = design_obj$rank,
    n_model_columns = design_obj$n_columns,
    residual_df = nrow(design) - design_obj$rank,
    condition_number = design_obj$condition_number,
    residual_disease_variance = design_obj$residual_disease_variance
  )]
  if (design_obj$rank != design_obj$n_columns ||
      !is.finite(design_obj$residual_disease_variance) ||
      design_obj$residual_disease_variance <= 1e-10 ||
      base_audit$residual_df < min_residual_df) {
    base_audit[, failure_reason := "rank_or_residual_information_failure"]
    return(list(
      effects = empty_effect_table(genes, cell_type, model_id,
                                   "rank_or_residual_information_failure"),
      audit = base_audit
    ))
  }

  ordered_counts <- counts[, as.character(meta$donor), drop = FALSE]
  dge <- DGEList(counts = ordered_counts)
  keep <- filterByExpr(dge, design = design)
  base_audit[, n_filter_by_expr := sum(keep)]
  if (sum(keep) < 100L) {
    base_audit[, failure_reason := "fewer_than_100_expressed_genes"]
    return(list(
      effects = empty_effect_table(genes, cell_type, model_id,
                                   "fewer_than_100_expressed_genes"),
      audit = base_audit
    ))
  }
  dge <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE], method = "TMM")
  v <- voom(dge, design, plot = FALSE)
  raw_fit <- lmFit(v, design)
  moderated <- eBayes(raw_fit, robust = TRUE)
  interval_fit <- treat(raw_fit, lfc = treat_lfc, robust = TRUE)
  k <- design_obj$coef_index
  tested_genes <- rownames(v$E)
  tested <- data.table(
    gene_symbol = tested_genes,
    logFC = moderated$coefficients[, k],
    moderated_se = moderated$stdev.unscaled[, k] * sqrt(moderated$s2.post),
    AveExpr = moderated$Amean,
    t = moderated$t[, k],
    pvalue = moderated$p.value[, k],
    B = moderated$lods[, k],
    treat_t = interval_fit$t[, k],
    treat_pvalue = interval_fit$p.value[, k]
  )
  tested[, qvalue := bh_adjust(pvalue)]
  tested[, treat_qvalue := bh_adjust(treat_pvalue)]

  ans <- data.table(gene_symbol = genes)
  ans[, filter_by_expr := gene_symbol %chin% tested_genes]
  ans <- merge(ans, tested, by = "gene_symbol", all.x = TRUE, sort = FALSE)
  ans[, `:=`(
    cell_type = cell_type, model_id = model_id, treat_lfc = treat_lfc,
    estimable = filter_by_expr & is.finite(logFC),
    failure_reason = fifelse(filter_by_expr, "", "not_filter_by_expr")
  )]
  setcolorder(ans, names(empty_effect_table(character(), cell_type, model_id, "")))
  base_audit[, `:=`(estimable = TRUE, failure_reason = "")]
  list(effects = ans, audit = base_audit)
}

message("[reference] loading approved metadata columns")
metadata_header <- names(fread(metadata_path, nrows = 0L))
required_metadata <- c(
  "sample", "dataset", "disease_stage_coarse", "exclude_stage_analysis",
  lineage_spec$cell_count_column
)
assert(all(required_metadata %in% metadata_header),
       "Extended metadata lacks required exact lineage counts: ",
       paste(setdiff(required_metadata, metadata_header), collapse = ";"))
run_meta <- fread(metadata_path, select = required_metadata)
run_meta[, exclude_stage_analysis := parse_flag(exclude_stage_analysis)]
for (column in lineage_spec$cell_count_column) {
  set(run_meta, j = column, value = as_numeric_zero(run_meta[[column]]))
}
donor_map <- build_srr_to_donor_map(project_root)
assert(length(donor_map) > 0L, "Authoritative donor map is empty")
run_meta[, donor := fifelse(sample %in% names(donor_map), donor_map[sample], sample)]

within_donor_audit <- run_meta[, .(
  n_source_records = .N,
  n_dataset_values = uniqueN(dataset[!is.na(dataset) & dataset != ""]),
  n_stage_values = uniqueN(disease_stage_coarse[
    !is.na(disease_stage_coarse) & disease_stage_coarse != ""
  ]),
  any_excluded = any(exclude_stage_analysis)
), by = donor]
bad_donors <- within_donor_audit[n_dataset_values > 1L | n_stage_values > 1L]
assert(nrow(bad_donors) == 0L,
       "Dataset/stage disagreement within biological donor: ",
       paste(bad_donors$donor, collapse = ";"))

donor_meta <- run_meta[, c(list(
  dataset = first_nonmissing(dataset),
  disease_stage_coarse = first_nonmissing(disease_stage_coarse),
  exclude_stage_analysis = any(exclude_stage_analysis),
  n_metadata_source_records = .N
), lapply(.SD, sum, na.rm = TRUE)), by = donor,
.SDcols = lineage_spec$cell_count_column]
donor_meta[, disease_binary := fcase(
  disease_stage_coarse == "Healthy", 0L,
  disease_stage_coarse %chin% c("Steatosis", "Steatohepatitis"), 1L,
  default = NA_integer_
)]

all_effects <- list()
all_audits <- list()
all_lodo <- list()
all_within <- list()
all_design_rows <- list()

for (lineage_index in seq_len(nrow(lineage_spec))) {
  cell_type <- lineage_spec$cell_type[[lineage_index]]
  count_column <- lineage_spec$cell_count_column[[lineage_index]]
  path <- file.path(pb_root, paste0(cell_type, "_pseudobulk.csv"))
  assert(file.exists(path), "Pseudobulk source absent: ", path)
  message("[reference] ", cell_type, " <- ", basename(path))
  raw <- fread(path, check.names = FALSE)
  genes <- as.character(raw[[1L]])
  counts <- as.matrix(raw[, -1L, with = FALSE])
  storage.mode(counts) <- "numeric"
  rownames(counts) <- genes
  if (anyDuplicated(genes)) counts <- rowsum(counts, group = genes, reorder = FALSE)
  counts <- collapse_counts_to_donor(counts, project_root)
  assert(!anyDuplicated(colnames(counts)), "Duplicate donor columns after collapse")
  libsize <- colSums(counts)

  design_rows <- donor_meta[, .(
    donor, dataset, disease_stage_coarse, disease_binary,
    exclude_stage_analysis, n_metadata_source_records,
    n_lineage_cells = get(count_column)
  )]
  design_rows[, `:=`(
    cell_type = cell_type,
    recipient_role = lineage_spec$recipient_role[[lineage_index]],
    count_matrix_present = donor %chin% colnames(counts),
    library_size = as.numeric(libsize[match(donor, names(libsize))])
  )]
  design_rows[, primary_eligible :=
    disease_stage_coarse %chin% primary_stages &
    !exclude_stage_analysis &
    dataset != source_uncertain_dataset &
    n_lineage_cells >= primary_min_cells &
    count_matrix_present & is.finite(library_size) & library_size > 0]
  design_rows[, min20_sensitivity_eligible :=
    disease_stage_coarse %chin% primary_stages &
    !exclude_stage_analysis &
    dataset != source_uncertain_dataset &
    n_lineage_cells >= sensitivity_min_cells &
    count_matrix_present & is.finite(library_size) & library_size > 0]
  design_rows[, include_source_uncertain_eligible :=
    disease_stage_coarse %chin% primary_stages &
    !exclude_stage_analysis &
    n_lineage_cells >= primary_min_cells &
    count_matrix_present & is.finite(library_size) & library_size > 0]
  design_rows[, primary_exclusion_reason := fcase(
    !disease_stage_coarse %chin% primary_stages, "stage_not_eligible",
    exclude_stage_analysis, "source_exclude_stage_analysis",
    dataset == source_uncertain_dataset, "source_independence_unconfirmed",
    n_lineage_cells < primary_min_cells, "fewer_than_50_lineage_cells",
    !count_matrix_present, "count_matrix_absent",
    !is.finite(library_size) | library_size <= 0, "zero_library",
    default = ""
  )]
  setorder(design_rows, cell_type, dataset, disease_stage_coarse, donor)
  all_design_rows[[cell_type]] <- design_rows

  model_specs <- list(
    primary = design_rows[primary_eligible == TRUE],
    min20_sensitivity = design_rows[min20_sensitivity_eligible == TRUE],
    include_source_uncertain_sensitivity = design_rows[
      include_source_uncertain_eligible == TRUE
    ]
  )
  for (model_id in names(model_specs)) {
    result <- fit_voom(counts, model_specs[[model_id]], cell_type, model_id)
    all_effects[[paste(cell_type, model_id, sep = "::")]] <- result$effects
    all_audits[[paste(cell_type, model_id, sep = "::")]] <- result$audit
  }

  primary_meta <- model_specs$primary
  for (held_out in sort(unique(primary_meta$dataset))) {
    lodo_meta <- primary_meta[dataset != held_out]
    lodo_id <- paste0("lodo_", held_out)
    result <- fit_voom(counts, lodo_meta, cell_type, lodo_id)
    all_lodo[[paste(cell_type, held_out, sep = "::")]] <- result$effects[, .(
      cell_type, held_out_dataset = held_out, gene_symbol,
      filter_by_expr, logFC, moderated_se, t, pvalue, qvalue,
      estimable, failure_reason
    )]
    audit <- copy(result$audit)
    audit[, model_id := lodo_id]
    all_audits[[paste(cell_type, lodo_id, sep = "::")]] <- audit
  }

  mixed <- primary_meta[, .(
    n_control = sum(disease_binary == 0L),
    n_disease = sum(disease_binary == 1L)
  ), by = dataset][n_control >= 2L & n_disease >= 2L]
  for (dataset_id in mixed$dataset) {
    cohort_meta <- primary_meta[dataset == dataset_id]
    result <- fit_voom(
      counts, cohort_meta, cell_type, paste0("within_", dataset_id),
      min_total = 4L, min_each = 2L, min_residual_df = 2L
    )
    all_within[[paste(cell_type, dataset_id, sep = "::")]] <- result$effects[, .(
      cell_type, dataset = dataset_id, gene_symbol,
      filter_by_expr, logFC, moderated_se, t, pvalue, qvalue,
      estimable, failure_reason
    )]
  }
  rm(raw, counts)
  gc(verbose = FALSE)
}

effects <- rbindlist(all_effects, use.names = TRUE, fill = TRUE)
audits <- rbindlist(all_audits, use.names = TRUE, fill = TRUE)
lodo <- rbindlist(all_lodo, use.names = TRUE, fill = TRUE)
within <- rbindlist(all_within, use.names = TRUE, fill = TRUE)
design_rows <- rbindlist(all_design_rows, use.names = TRUE, fill = TRUE)

primary <- effects[model_id == "primary"]
lodo_summary <- lodo[estimable & is.finite(logFC), .(
  n_lodo_estimable = .N,
  n_lodo_positive = sum(logFC > 0),
  n_lodo_negative = sum(logFC < 0)
), by = .(cell_type, gene_symbol)]
within_summary <- within[estimable & is.finite(logFC), .(
  n_within_dataset_estimable = .N,
  n_within_positive = sum(logFC > 0),
  n_within_negative = sum(logFC < 0)
), by = .(cell_type, gene_symbol)]
reference <- merge(primary, lodo_summary, by = c("cell_type", "gene_symbol"), all.x = TRUE)
reference <- merge(reference, within_summary,
                   by = c("cell_type", "gene_symbol"), all.x = TRUE)
for (column in c(
  "n_lodo_estimable", "n_lodo_positive", "n_lodo_negative",
  "n_within_dataset_estimable", "n_within_positive", "n_within_negative"
)) {
  set(reference, which(is.na(reference[[column]])), column, 0L)
}
reference[, static_mitochondrial := grepl("^MT-", gene_symbol)]
reference[, static_ribosomal := grepl("^(RPL|RPS)[0-9A-Z-]*$", gene_symbol)]
reference[, static_annotation_identity := gene_symbol %chin% identity_markers]
reference[, static_nuisance :=
  static_mitochondrial | static_ribosomal | static_annotation_identity]
reference[, primary_direction := sign(logFC)]
reference[, lodo_direction_agreement_fraction := fifelse(
  n_lodo_estimable > 0 & primary_direction != 0,
  fifelse(primary_direction > 0, n_lodo_positive, n_lodo_negative) /
    n_lodo_estimable,
  NA_real_
)]
reference[, within_direction_agreement_fraction := fifelse(
  n_within_dataset_estimable > 0 & primary_direction != 0,
  fifelse(primary_direction > 0, n_within_positive, n_within_negative) /
    n_within_dataset_estimable,
  NA_real_
)]
reference[, primary_reference_eligible :=
  estimable & is.finite(logFC) & logFC != 0 & !static_nuisance &
  n_lodo_estimable >= 3L & lodo_direction_agreement_fraction >= 0.75]
reference[, strict_treat_reference_eligible :=
  primary_reference_eligible & is.finite(treat_qvalue) & treat_qvalue < 0.05]
reference[, within_dataset_direction_sensitivity_pass :=
  n_within_dataset_estimable >= 2L & within_direction_agreement_fraction >= 1]
reference[, primary_raw_loading := fifelse(primary_reference_eligible, logFC, 0)]
reference[, sign_only_raw_loading := fifelse(
  primary_reference_eligible, as.numeric(sign(logFC)), 0
)]
reference[, strict_treat_raw_loading := fifelse(
  strict_treat_reference_eligible, logFC, 0
)]
reference[, dynamic_exclusions_required := paste(
  "remove the perturbed cis target and every program containing it;",
  "remove the preregistered guide-response gene set before experimental outcome scoring"
)]
setorder(reference, cell_type, gene_symbol)

gate <- audits[model_id == "primary", .(
  cell_type, model_id, n_donors, n_controls, n_disease, n_datasets,
  n_mixed_datasets, design_rank, n_model_columns, residual_df,
  condition_number, residual_disease_variance, n_filter_by_expr,
  estimable, failure_reason
)]
gate[, `:=`(
  required_min_controls = 5L,
  required_min_disease = 5L,
  required_mixed_datasets = 2L,
  source_gate_pass = estimable & n_controls >= 5L & n_disease >= 5L &
    n_mixed_datasets >= 2L & design_rank == n_model_columns & residual_df >= 5L
)]
gate[!source_gate_pass & failure_reason == "", failure_reason :=
  "fewer_than_two_independent_within_dataset_contrasts"]
assert(all(gate$source_gate_pass),
       "At least one primary lineage failed the frozen source gate: ",
       paste(gate[!source_gate_pass, cell_type], collapse = ";"))

nuisance <- unique(reference[
  static_nuisance == TRUE,
  .(gene_symbol, static_mitochondrial, static_ribosomal,
    static_annotation_identity)
])
setorder(nuisance, gene_symbol)

dir.create(candidate_parent, recursive = TRUE, showWarnings = FALSE)
temp_root <- tempfile(pattern = paste0(".", candidate_id, "."), tmpdir = candidate_parent)
dir.create(temp_root)
on.exit(if (dir.exists(temp_root)) unlink(temp_root, recursive = TRUE), add = TRUE)

output_objects <- list(
  donor_lineage_design = design_rows,
  lineage_disease_coefficients = effects,
  lineage_design_audit = audits,
  lineage_disease_lodo = lodo,
  lineage_within_dataset_effects = within,
  frozen_lineage_disease_reference = reference,
  static_nuisance_gene_registry = nuisance,
  source_gate_status = gate,
  donor_metadata_collapse_audit = within_donor_audit
)
output_paths <- list()
for (name in names(output_objects)) {
  path <- file.path(temp_root, paste0(name, ".tsv"))
  atomic_fwrite(output_objects[[name]], path)
  output_paths[[name]] <- path
}

pairing_paths <- file.path(project_root, .DONOR_PAIRING_FILES)
pb_paths <- file.path(pb_root, paste0(lineage_spec$cell_type, "_pseudobulk.csv"))
input_manifest <- rbindlist(list(
  data.table(
    source_id = paste0("pseudobulk_", lineage_spec$cell_type),
    source_path = sub(paste0("^", project_root, "/"), "", pb_paths),
    size_bytes = file.info(pb_paths)$size,
    sha256 = vapply(pb_paths, sha256, character(1L)),
    role = "run-keyed raw-count lineage pseudobulk"
  ),
  data.table(
    source_id = "run_metadata_with_exact_lineage_cell_counts",
    source_path = sub(paste0("^", project_root, "/"), "", metadata_path),
    size_bytes = file.info(metadata_path)$size,
    sha256 = sha256(metadata_path),
    role = "approved stage fields and exact per-run lineage cell counts"
  ),
  data.table(
    source_id = "donor_collapse_implementation",
    source_path = sub(paste0("^", project_root, "/"), "", collapse_path),
    size_bytes = file.info(collapse_path)$size,
    sha256 = sha256(collapse_path),
    role = "authoritative run-to-donor collapse implementation"
  ),
  data.table(
    source_id = paste0("donor_pairing_", names(.DONOR_PAIRING_FILES)),
    source_path = .DONOR_PAIRING_FILES,
    size_bytes = file.info(pairing_paths)$size,
    sha256 = vapply(pairing_paths, sha256, character(1L)),
    role = "authoritative run-to-biological-donor source"
  ),
  data.table(
    source_id = "original_annotation_marker_source",
    source_path = sub(paste0("^", project_root, "/"), "", annotation_path),
    size_bytes = file.info(annotation_path)$size,
    sha256 = sha256(annotation_path),
    role = "version-controlled identity-marker exclusions"
  ),
  data.table(
    source_id = "sealed_baseline_lineage_reference",
    source_path = sub(paste0("^", project_root, "/"), "", lineage_reference_seal),
    size_bytes = file.info(lineage_reference_seal)$size,
    sha256 = sha256(lineage_reference_seal),
    role = "outcome-blind upstream lineage observability seal"
  )
), use.names = TRUE)
manifest_path <- file.path(temp_root, "lineage_disease_input_manifest.tsv")
atomic_fwrite(input_manifest, manifest_path)
output_paths$lineage_disease_input_manifest <- manifest_path

seal <- list(
  status = "sealed_lineage_resolved_disease_reference",
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  candidate_id = candidate_id,
  biological_unit = "biological donor",
  primary_contrast = "Steatosis_or_Steatohepatitis_minus_Healthy",
  model = "lineage raw-count pseudobulk; TMM; limma-voom; dataset fixed effect",
  primary_min_lineage_cells_per_donor = primary_min_cells,
  sensitivity_min_lineage_cells_per_donor = sensitivity_min_cells,
  primary_excluded_dataset = source_uncertain_dataset,
  treat_interval_log2fc = treat_lfc,
  primary_lineages = lineage_spec[recipient_role == "primary_organoid", cell_type],
  optional_lineages = lineage_spec[recipient_role == "optional_recipient", cell_type],
  reference_rule = paste(
    "finite primary coefficient; static nuisance excluded; at least three",
    "estimable leave-one-dataset-out fits; at least 75% direction agreement"
  ),
  primary_score_contract = paste(
    "within-lineage experimental pseudobulk expression standardized to",
    "protective/control samples; project onto primary_raw_loading and normalize",
    "by retained absolute loading after target and preregistered guide-response exclusions"
  ),
  failed_plan43_or_plan44_score_reused = FALSE,
  inferred_f_stage_used = FALSE,
  run_level_rows_used_as_replicates = FALSE,
  experimental_outcomes_inspected = FALSE,
  experimental_targets_frozen = FALSE,
  static_identity_markers = identity_markers,
  dynamic_exclusions_pending = c("perturbed_cis_target", "preregistered_guide_response_genes"),
  output_sha256 = lapply(output_paths, sha256)
)
seal_path <- file.path(temp_root, "LINEAGE_DISEASE_REFERENCE_SEALED.json")
writeLines(toJSON(seal, auto_unbox = TRUE, pretty = TRUE), seal_path)

if (!file.rename(temp_root, candidate_root)) stop("Atomic candidate promotion failed")
message(
  "[reference] sealed ", nrow(reference), " lineage-gene rows; ",
  sum(reference$primary_reference_eligible), " primary loading rows; targets not frozen"
)
