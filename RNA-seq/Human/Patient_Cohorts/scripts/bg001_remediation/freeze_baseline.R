#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: freeze_baseline.R PROJECT_ROOT RUN_ROOT", call. = FALSE)
project_root <- normalizePath(args[[1]], mustWork = TRUE)
run_root <- normalizePath(args[[2]], mustWork = TRUE)
if (!file.exists(file.path(run_root, ".bg001_candidate_root"))) stop("Missing candidate sentinel", call. = FALSE)
snapshot_root <- normalizePath(file.path(run_root, "source_snapshot"), mustWork = TRUE)

patient_root <- file.path(project_root, "RNA-seq/Human/Patient_Cohorts")
source_rdir <- file.path(patient_root, "analysis/integration/results/integration")
source_qc <- file.path(patient_root, "analysis/integration/qc/sample_qc_report.csv")
frozen <- file.path(run_root, "frozen_sets")
dir.create(frozen, recursive = TRUE, showWarnings = FALSE)

# Refuse to freeze any live baseline object unless the run snapshot contains the
# exact 78-panel rendered-label sidecars and every frozen PDF/text binding still
# validates. This executes the snapshot's validator, not concurrent worktree code.
rendered_validator <- file.path(
  snapshot_root,
  "RNA-seq/Human/Patient_Cohorts/scripts/bg001_remediation/validate_rendered_figure_labels.py"
)
if (!file.exists(rendered_validator)) stop("Missing frozen rendered-label validator", call. = FALSE)
rendered_check <- suppressWarnings(system2(
  "python3",
  c(shQuote(rendered_validator), "--project-root", shQuote(snapshot_root)),
  stdout = TRUE,
  stderr = TRUE
))
rendered_status <- attr(rendered_check, "status")
if (is.null(rendered_status)) rendered_status <- 0L
writeLines(
  c(sprintf("exit_status=%d", rendered_status), rendered_check),
  file.path(frozen, "rendered_figure_label_validation.txt")
)
if (rendered_status != 0L) {
  stop("Frozen rendered figure-label validation failed before baseline freeze", call. = FALSE)
}

file_sha256 <- function(path) {
  output <- system2("sha256sum", shQuote(path), stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) stop("sha256sum failed for ", path, call. = FALSE)
  strsplit(output[[1L]], "[[:space:]]+")[[1L]][[1L]]
}

file_identity <- function(path) {
  output <- system2(
    "stat", c("-c", shQuote("%d\t%i\t%s\t%Y\t%f"), "--", shQuote(path)),
    stdout = TRUE, stderr = TRUE
  )
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) stop("stat failed for ", path, call. = FALSE)
  output[[1L]]
}

copy_exact <- function(source, destination) {
  if (!file.exists(source) || nzchar(Sys.readlink(source)) || isTRUE(file.info(source)$isdir)) {
    stop("Baseline source is missing, symlinked, or non-regular: ", source, call. = FALSE)
  }
  # Sys.readlink() returns NA for a path that does not exist, and nzchar(NA) is
  # TRUE by default, so a bare nzchar() here fired on EVERY absent destination —
  # i.e. it refused exactly the normal case. keepNA = TRUE + isTRUE() makes the
  # symlink test NA-safe while preserving the intent.
  if (file.exists(destination) ||
      isTRUE(nzchar(Sys.readlink(destination), keepNA = TRUE))) {
    stop("Refusing existing frozen destination: ", destination, call. = FALSE)
  }
  source_identity_before <- file_identity(source)
  source_hash_before <- file_sha256(source)
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(source, destination, overwrite = FALSE, copy.date = TRUE)) {
    stop("Failed to freeze ", source, call. = FALSE)
  }
  source_identity_after <- file_identity(source)
  source_hash_after <- file_sha256(source)
  if (!identical(source_identity_before, source_identity_after) ||
      !identical(source_hash_before, source_hash_after)) {
    stop("Baseline source changed while it was being frozen: ", source, call. = FALSE)
  }
  if (!file.exists(destination) || nzchar(Sys.readlink(destination)) ||
      isTRUE(file.info(destination)$isdir) || file_sha256(destination) != source_hash_before) {
    stop("Frozen destination is non-regular or differs from its source: ", destination, call. = FALSE)
  }
  destination
}

# Freeze the exact live objects before deriving any baseline sets. All analysis
# arms read these copies, never the live canonical paths.
locked_qc_path <- copy_exact(source_qc, file.path(frozen, "locked_sample_qc_report.csv"))
locked_meta_path <- copy_exact(
  file.path(source_rdir, "meta_matched.rds"),
  file.path(frozen, "locked_meta_matched.rds")
)
locked_dge_path <- copy_exact(
  file.path(source_rdir, "merged_dge.rds"),
  file.path(frozen, "locked_merged_dge.rds")
)
canonical_deg_path <- copy_exact(
  file.path(source_rdir, "canonical_deg_results.csv"),
  file.path(frozen, "canonical_deg_results.csv")
)
evidence_source <- file.path(snapshot_root, "RNA-seq/results/manuscript_release/2026-07-15-r2/evidence_class_table.tsv")
evidence_frozen <- copy_exact(evidence_source, file.path(frozen, "evidence_class_table_2026-07-15-r2.tsv"))

dge <- readRDS(locked_dge_path)
meta <- as.data.table(readRDS(locked_meta_path))
qc <- fread(locked_qc_path)
deg <- fread(canonical_deg_path)
cfg <- yaml::read_yaml(file.path(snapshot_root, "config/human_datasets.yaml"))$datasets
mega <- names(Filter(function(x) isTRUE(x$de$include_in_mega), cfg))
keep_mega <- dge$samples$dataset %in% mega
canonical_samples <- colnames(dge)[keep_mega]
if (length(canonical_samples) != 846L || nrow(dge) != 27638L) stop("Frozen canonical dimensions changed", call. = FALSE)
if (nrow(qc) != ncol(readRDS(file.path(source_rdir, "merged_counts_raw.rds")))) {
  stop("Canonical QC and raw-count sample dimensions disagree", call. = FALSE)
}
substrate_census <- data.table(
  matched_samples = nrow(qc),
  pass_technical_samples = sum(qc$pass_technical %in% TRUE),
  dge_samples = ncol(dge),
  canonical_model_samples = length(canonical_samples)
)
if (substrate_census$matched_samples != 1281L ||
    substrate_census$pass_technical_samples != 1260L ||
    substrate_census$dge_samples != 1260L ||
    substrate_census$canonical_model_samples != 846L) {
  stop("Authoritative 2026-08-06 integration-substrate census changed", call. = FALSE)
}
fwrite(substrate_census, file.path(frozen, "integration_substrate_census.tsv"), sep = "\t")
treat <- deg[is.finite(treat_fdr) & treat_fdr < 0.05]
if (nrow(treat) != 1918L || sum(treat$logFC > 0) != 1419L || sum(treat$logFC < 0) != 499L) {
  stop("Frozen canonical TREAT contract changed", call. = FALSE)
}

# Freeze every current read-mode matrix used by R0 and the unaffected arms.
read_count_rows <- list()
for (dataset in names(cfg)) {
  source_count <- file.path(patient_root, cfg[[dataset]]$counts_path)
  frozen_count <- file.path(frozen, "read_counts", dataset, "gene_counts.txt")
  copy_exact(source_count, frozen_count)
  source_summary <- paste0(source_count, ".summary")
  frozen_summary <- paste0(frozen_count, ".summary")
  if (file.exists(source_summary)) copy_exact(source_summary, frozen_summary)
  read_count_rows[[dataset]] <- data.table(
    dataset = dataset,
    source_count_path = normalizePath(source_count, mustWork = TRUE),
    frozen_count_path = normalizePath(frozen_count, mustWork = TRUE),
    source_summary_path = if (file.exists(source_summary)) normalizePath(source_summary, mustWork = TRUE) else "",
    frozen_summary_path = if (file.exists(frozen_summary)) normalizePath(frozen_summary, mustWork = TRUE) else ""
  )
}
read_count_manifest <- rbindlist(read_count_rows)
fwrite(read_count_manifest, file.path(run_root, "manifests/read_count_sources.tsv"), sep = "\t")
fwrite(
  read_count_manifest[, .(dataset, count_path = frozen_count_path)],
  file.path(run_root, "manifests/read_count_overrides.tsv"),
  sep = "\t"
)

# Exact sample, gene, QC, sex, design, and normalization baselines.
fwrite(qc[pass_technical == TRUE, .(sample_id, dataset, pass_technical)], file.path(frozen, "locked_pass_technical.tsv"), sep = "\t")
fwrite(meta[, .(sample_id, dataset, group_binary, inferred_sex, sex_final)], file.path(frozen, "locked_metadata.tsv"), sep = "\t")
fwrite(data.table(ordinal = seq_along(canonical_samples), sample_id = canonical_samples), file.path(frozen, "canonical_samples.tsv"), sep = "\t")
fwrite(data.table(ordinal = seq_len(nrow(dge)), gene = rownames(dge)), file.path(frozen, "canonical_genes.tsv"), sep = "\t")
fwrite(treat[, .(gene, symbol, logFC, treat_fdr)], file.path(frozen, "canonical_treat_genes.tsv"), sep = "\t")

model_meta <- meta[match(canonical_samples, sample_id)]
if (anyNA(model_meta$sample_id) || anyNA(model_meta$inferred_sex)) stop("Cannot freeze baseline model metadata", call. = FALSE)
model_info <- data.frame(
  group_binary = factor(dge$samples$group_binary[keep_mega], levels = c("Control", "Disease")),
  dataset = factor(dge$samples$dataset[keep_mega]),
  inferred_sex = factor(model_meta$inferred_sex)
)
baseline_design <- model.matrix(~ dataset + inferred_sex + group_binary, data = model_info)
if (qr(baseline_design)$rank != ncol(baseline_design) || !"group_binaryDisease" %in% colnames(baseline_design)) {
  stop("Frozen baseline design is not full rank with an estimable disease coefficient", call. = FALSE)
}
fwrite(as.data.table(baseline_design, keep.rownames = "sample_id")[, sample_id := canonical_samples],
       file.path(frozen, "canonical_model_design.tsv"), sep = "\t")
baseline_norm <- as.data.table(dge$samples, keep.rownames = "sample_id")
baseline_norm[, effective_library_size := lib.size * norm.factors]
fwrite(baseline_norm, file.path(frozen, "canonical_normalization_factors.tsv"), sep = "\t")

# Finite-SuSiE comparison set (held fixed against concurrent BG-031 work).
evidence <- fread(evidence_frozen)
susie <- unique(evidence[in_susie_primary == TRUE & bulk_tested == TRUE & !is.na(ensembl_bulk), .(
  gene = ensembl_bulk,
  symbol,
  max_susie_pp4_all
)], by = "gene")
if (nrow(susie) != 447L) stop("Frozen finite-SuSiE primary set is not 447 genes", call. = FALSE)
fwrite(susie, file.path(frozen, "finite_susie_447.tsv"), sep = "\t")

# Freeze and scan the active manuscript/figure/benchmark/target sources. Only
# genes testable in the old 27,638-gene universe are protected; unmatched names
# are retained separately so the scope is auditable rather than silently lost.
text_sources <- unique(c(
  file.path(snapshot_root, "docs/archive/documentation_consolidation_2026-08-11/originals/docs/manuscript/NUMBERS.md"),
  file.path(snapshot_root, "docs/archive/documentation_consolidation_2026-08-11/originals/docs/paper_outline.md"),
  list.files(file.path(snapshot_root, "docs/archive/documentation_consolidation_2026-08-11/originals/docs/manuscript/working"), pattern = "[.]md$", full.names = TRUE, recursive = TRUE),
  list.files(file.path(snapshot_root, "scripts/figures"), pattern = "[.]R$", full.names = TRUE, recursive = TRUE),
  Sys.glob(file.path(snapshot_root, "RNA-seq/46d*.R"))
))
table_sources <- c(
  positive_control = file.path(snapshot_root, "RNA-seq/results/validation/positive_control_validation.csv"),
  triple_convergence = file.path(snapshot_root, "RNA-seq/results/convergence/triple_convergence_targets.csv"),
  inhibitor_targets = file.path(snapshot_root, "RNA-seq/results/multi_evidence/convergence_evidence_inhibitor_target_candidates.csv"),
  heldout_govaere = file.path(snapshot_root, "data/published_gene_panels/govaere_2020_panel.tsv"),
  heldout_niddk = file.path(snapshot_root, "data/published_gene_panels/niddk_pipeline_2024.tsv"),
  heldout_opentargets = file.path(snapshot_root, "data/published_gene_panels/opentargets_masld_2025.tsv"),
  panel4c_fixed = file.path(snapshot_root, "Analysis/Multimodal_Program_Projection/config/panel4c_fixed_rows.tsv"),
  curated_named = file.path(snapshot_root, "RNA-seq/Human/Patient_Cohorts/scripts/bg001_remediation/protected_named_genes.tsv"),
  curated_figure = file.path(snapshot_root, "RNA-seq/Human/Patient_Cohorts/scripts/bg001_remediation/protected_figure_genes.tsv"),
  rendered_figure_scope = file.path(snapshot_root, "RNA-seq/Human/Patient_Cohorts/scripts/bg001_remediation/active_figure_panel_scope.tsv"),
  rendered_figure_roster = file.path(snapshot_root, "RNA-seq/Human/Patient_Cohorts/scripts/bg001_remediation/rendered_figure_label_roster.tsv"),
  rendered_figure_inventory = file.path(snapshot_root, "RNA-seq/Human/Patient_Cohorts/scripts/bg001_remediation/active_figure_panel_inventory.tsv"),
  rendered_figure_tokens = file.path(snapshot_root, "RNA-seq/Human/Patient_Cohorts/scripts/bg001_remediation/rendered_figure_label_tokens.tsv")
)
all_sources <- unique(c(text_sources[file.exists(text_sources)], unname(table_sources[file.exists(table_sources)])))
source_rows <- list()
frozen_source_paths <- character()
for (source in all_sources) {
  relative <- substring(normalizePath(source, mustWork = TRUE), nchar(snapshot_root) + 2L)
  destination <- file.path(frozen, "protected_sources", relative)
  copy_exact(source, destination)
  source_rows[[length(source_rows) + 1L]] <- data.table(
    source_path = file.path(project_root, relative),
    snapshot_path = source,
    frozen_path = destination
  )
  frozen_source_paths[[source]] <- destination
}
fwrite(rbindlist(source_rows), file.path(run_root, "manifests/protected_source_manifest.tsv"), sep = "\t")

gene_meta <- fread(file.path(snapshot_root, "data/gencode_v49_gene_metadata.tsv.gz"))
current_map <- unique(data.table(
  gene = rownames(dge),
  ensembl_base = sub("[.][0-9]+$", "", rownames(dge))
))
gene_meta <- unique(gene_meta[, .(gene = gene_id, symbol = gene_name, ensembl_base)], by = "gene")
current_map <- merge(current_map, gene_meta[, .(ensembl_base, symbol)], by = "ensembl_base", all.x = TRUE)
setcolorder(current_map, c("gene", "ensembl_base", "symbol"))

validate_curated_inventory <- function(source) {
  tab <- fread(frozen_source_paths[[source]])
  expected_columns <- c("symbol", "source_path", "source_line", "context", "baseline_status")
  if (!identical(names(tab), expected_columns) || anyDuplicated(tab$symbol) ||
      anyNA(tab[, ..expected_columns]) || any(!nzchar(tab$symbol))) {
    stop("Malformed curated protected-gene inventory: ", source, call. = FALSE)
  }
  official <- unique(gene_meta$symbol)
  if (any(!tab$symbol %in% official)) stop("Curated inventory contains a non-GENCODE-v49 symbol: ", source, call. = FALSE)
  expected_status <- ifelse(tab$symbol %in% current_map$symbol, "testable", "NOT_TESTABLE")
  if (!identical(tab$baseline_status, expected_status)) {
    stop("Curated inventory baseline status differs from the frozen DGE universe: ", source, call. = FALSE)
  }
  for (index in seq_len(nrow(tab))) {
    relative <- tab$source_path[[index]]
    if (grepl("^/|(^|/)[.][.](/|$)", relative)) stop("Curated source path escapes snapshot: ", relative, call. = FALSE)
    path <- file.path(snapshot_root, relative)
    if (!file.exists(path)) stop("Curated source is absent from the snapshot: ", relative, call. = FALSE)
    source_lines <- readLines(path, warn = FALSE)
    line_number <- as.integer(tab$source_line[[index]])
    if (!is.finite(line_number) || line_number < 1L || line_number > length(source_lines) ||
        !grepl(tab$symbol[[index]], source_lines[[line_number]], fixed = TRUE)) {
      stop("Curated source line does not contain its symbol: ", relative, ":", line_number, call. = FALSE)
    }
  }
  invisible(tab)
}
validate_curated_inventory(table_sources[["curated_named"]])
validate_curated_inventory(table_sources[["curated_figure"]])

protected_parts <- list(
  treat[, .(gene, symbol, protected_reason = "canonical_treat", protected_source = "canonical_deg_results.csv")],
  current_map[ensembl_base %in% sub("[.][0-9]+$", "", susie$gene),
              .(gene, symbol, protected_reason = "finite_susie_primary", protected_source = "evidence_class_table_2026-07-15-r2.tsv")]
)
named_candidates <- list()

add_symbols <- function(symbols, reason, source) {
  symbols <- unique(symbols[!is.na(symbols) & nzchar(symbols)])
  if (!length(symbols)) return(invisible(NULL))
  named_candidates[[length(named_candidates) + 1L]] <<- data.table(symbol = symbols, protected_reason = reason, protected_source = source)
  hits <- current_map[!is.na(symbol) & symbol %in% symbols]
  if (nrow(hits)) {
    protected_parts[[length(protected_parts) + 1L]] <<- hits[, .(
      gene, symbol, protected_reason = reason, protected_source = source
    )]
  }
}

read_heldout_panel <- function(path) {
  lines <- readLines(path, warn = FALSE)
  header_lines <- which(startsWith(lines, "gene_symbol\t") | lines == "gene_symbol")
  if (length(header_lines) != 1L) {
    stop("Expected exactly one gene_symbol table header in held-out panel: ", path, call. = FALSE)
  }
  tab <- fread(path, skip = header_lines[[1L]] - 1L)
  if (!"gene_symbol" %in% names(tab) || !nrow(tab) ||
      anyNA(tab$gene_symbol) || any(!nzchar(tab$gene_symbol))) {
    stop("Malformed or empty gene_symbol column in held-out panel: ", path, call. = FALSE)
  }
  tab
}

for (source_name in names(table_sources)) {
  source <- table_sources[[source_name]]
  if (!file.exists(source)) next
  tab <- if (grepl("^heldout_", source_name)) {
    read_heldout_panel(frozen_source_paths[[source]])
  } else {
    fread(frozen_source_paths[[source]])
  }
  symbol_columns <- intersect(c("gene", "symbol", "human_symbol", "gene_name", "gene_symbol", "tf_name"), names(tab))
  for (column in symbol_columns) add_symbols(tab[[column]], paste0("active_", source_name), basename(source))
}

# Dynamic rendered labels are the only PDF-derived rows added here. Static
# labels already have exact source-line proofs in protected_figure_genes.tsv;
# collision rows are deliberately excluded. add_symbols() retains official
# names absent from the old DGE universe in named_genes_not_baseline_testable.tsv.
rendered_scope <- fread(frozen_source_paths[[table_sources[["rendered_figure_scope"]]]])
rendered_tokens <- fread(frozen_source_paths[[table_sources[["rendered_figure_tokens"]]]])
if (!identical(names(rendered_scope), c("panel_id", "pdf_path", "label_mode", "classification_reason")) ||
    !identical(names(rendered_tokens), c("panel_id", "pdf_path", "token", "classification", "reason", "proof_source_path", "proof_source_line"))) {
  stop("Malformed frozen rendered figure-label sidecars", call. = FALSE)
}
dynamic_panels <- rendered_scope[label_mode == "dynamic_rendered", panel_id]
if (length(dynamic_panels) != 24L || anyDuplicated(dynamic_panels)) {
  stop("Frozen rendered-label scope is not the reviewed 24-panel dynamic roster", call. = FALSE)
}
dynamic_gene_rows <- rendered_tokens[
  panel_id %in% dynamic_panels & classification == "gene_label"
]
if (!nrow(dynamic_gene_rows) || any(!dynamic_gene_rows$pdf_path %in% rendered_scope$pdf_path)) {
  stop("Frozen dynamic rendered-label gene rows are empty or unscoped", call. = FALSE)
}
add_symbols(
  dynamic_gene_rows$token,
  "active_rendered_figure_label",
  basename(table_sources[["rendered_figure_tokens"]])
)

protected_long <- rbindlist(protected_parts, use.names = TRUE, fill = TRUE)
protected <- protected_long[, .(
  symbol = if (all(is.na(symbol))) NA_character_ else symbol[which(!is.na(symbol))[1L]],
  protected_reason = paste(sort(unique(protected_reason)), collapse = "|"),
  protected_source = paste(sort(unique(protected_source)), collapse = "|")
), by = gene]
fwrite(protected, file.path(frozen, "protected_genes.tsv"), sep = "\t")
candidate_names <- unique(rbindlist(named_candidates, use.names = TRUE, fill = TRUE))
unmatched <- candidate_names[!symbol %in% current_map$symbol]
fwrite(unmatched, file.path(frozen, "named_genes_not_baseline_testable.tsv"), sep = "\t")

writeLines(capture.output(sessionInfo()), file.path(frozen, "baseline_session_info.txt"))
cat(sprintf(
  "Frozen baseline: 846 canonical samples, 27,638 genes, 1,918 TREAT genes, 447 finite-SuSiE genes, %d protected baseline-testable genes\n",
  nrow(protected)
))
