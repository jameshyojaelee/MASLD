#!/usr/bin/env Rscript
# Intersect pair-specific shared signals with exact source-lineage adult-liver
# accessibility. This dependency-gated stage still cannot freeze a target.

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(GenomicRanges)
  library(rtracklayer)
})

project_root <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
candidate_id <- Sys.getenv("PLAN45_CANDIDATE_ID", unset = "")
if (!grepl("^source-independent-risk-state-relay-[A-Za-z0-9._-]+$", candidate_id)) {
  stop("PLAN45_CANDIDATE_ID must name a new routing candidate")
}
candidate_root <- file.path(
  project_root,
  "Analysis/Multimodal_Program_Projection/candidates",
  candidate_id
)
if (dir.exists(candidate_root)) stop("Refusing to overwrite candidate: ", candidate_root)

resolve_candidate <- function(variable) {
  raw <- Sys.getenv(variable, unset = "")
  if (!nzchar(raw)) stop(variable, " is required")
  path <- if (grepl("^/", raw)) raw else file.path(project_root, raw)
  normalized <- normalizePath(path, mustWork = TRUE)
  allowed <- normalizePath(
    file.path(project_root, "Analysis/Multimodal_Program_Projection/candidates"),
    mustWork = TRUE
  )
  if (!startsWith(normalized, paste0(allowed, "/"))) {
    stop(variable, " escapes the candidate root")
  }
  normalized
}

orientation_root <- resolve_candidate("PLAN45_ORIENTATION_ROOT")
lineage_root <- resolve_candidate("PLAN45_LINEAGE_ROOT")
accessibility_root <- resolve_candidate("PLAN45_ACCESSIBILITY_ROOT")

sha256 <- function(path) {
  output <- system2("sha256sum", args = path, stdout = TRUE, stderr = TRUE)
  if (length(output) != 1L) stop("sha256sum failed: ", path)
  strsplit(output, "[[:space:]]+")[[1L]][[1L]]
}

as_flag <- function(x) tolower(as.character(x)) %chin% c("true", "t", "1", "yes")

atomic_fwrite <- function(x, path) {
  tmp <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  fwrite(x, tmp, sep = "\t", quote = FALSE, na = "")
  if (!file.rename(tmp, path)) stop("Atomic rename failed: ", path)
}

require_hash <- function(path, expected) {
  if (!file.exists(path) || file.info(path)$size <= 0L) stop("Missing input: ", path)
  if (!identical(sha256(path), expected)) stop("Input hash mismatch: ", path)
}

# Validate sealed upstream products before loading their tables.
orientation_seal_path <- file.path(orientation_root, "ORIENTATION_RELEASE.json")
lineage_seal_path <- file.path(lineage_root, "LINEAGE_REFERENCE_SEALED.json")
access_seal_path <- file.path(accessibility_root, "ACCESSIBILITY_REFERENCE_SEALED.json")
for (path in c(orientation_seal_path, lineage_seal_path, access_seal_path)) {
  if (!file.exists(path)) stop("Missing upstream seal: ", path)
}
orientation_seal <- fromJSON(orientation_seal_path, simplifyVector = FALSE)
lineage_seal <- fromJSON(lineage_seal_path, simplifyVector = FALSE)
access_seal <- fromJSON(access_seal_path, simplifyVector = FALSE)
if (!identical(orientation_seal$status, "orientation_evidence_frozen_targets_not_selected")) {
  stop("Orientation release status is not eligible")
}
if (!identical(lineage_seal$status, "donor_collapsed_lineage_observability_reference")) {
  stop("Lineage release status is not eligible")
}
if (!identical(access_seal$status, "accessibility_routing_reference_targets_not_selected")) {
  stop("Accessibility release status is not eligible")
}
if (isTRUE(orientation_seal$experimental_targets_frozen) ||
    isTRUE(lineage_seal$experimental_targets_frozen) ||
    isTRUE(access_seal$experimental_targets_frozen)) {
  stop("An upstream reference improperly froze targets")
}

orientation_path <- file.path(orientation_root, "shared_signal_orientation.tsv")
worklist_path <- file.path(orientation_root, "orientation_worklist.tsv")
pair_manifest_path <- file.path(orientation_root, "pair_export_manifest.tsv")
lineage_path <- file.path(lineage_root, "donor_lineage_gene_summary.tsv")
access_manifest_path <- file.path(accessibility_root, "accessibility_source_manifest.tsv")
access_map_path <- file.path(accessibility_root, "lineage_accessibility_map.tsv")
require_hash(orientation_path, orientation_seal$shared_signal_orientation_sha256)
require_hash(pair_manifest_path, orientation_seal$pair_export_manifest_sha256)
require_hash(
  lineage_path,
  lineage_seal$output_sha256$donor_lineage_gene_summary
)
require_hash(
  access_manifest_path,
  access_seal$output_sha256$accessibility_source_manifest.tsv
)
require_hash(
  access_map_path,
  access_seal$output_sha256$lineage_accessibility_map.tsv
)

orientation <- fread(orientation_path)
worklist <- fread(worklist_path)
pair_manifest <- fread(pair_manifest_path)
lineage <- fread(lineage_path)
access_manifest <- fread(access_manifest_path)
access_map <- fread(access_map_path)
if (anyDuplicated(orientation$orientation_uid) || anyDuplicated(worklist$orientation_uid)) {
  stop("Duplicate orientation UID")
}
if (!setequal(orientation$orientation_uid, worklist$orientation_uid)) {
  stop("Orientation/worklist universe mismatch")
}
if (anyDuplicated(lineage$gene_symbol)) stop("Duplicate lineage gene summary rows")
if (nrow(access_map) != 16L || sum(access_map$mapping_status == "exact") != 9L) {
  stop("Accessibility lineage-map drift")
}

# Verify the pair-export variant payload hashes, then assemble all prespecified
# shared-95 variants. Failed upstream orientations remain represented in the
# locus table but do not enter the accessibility gate.
variant_manifest <- pair_manifest[role == "variant_posterior"]
if (nrow(variant_manifest) != nrow(worklist)) stop("Variant manifest is incomplete")
variant_parts <- vector("list", nrow(variant_manifest))
for (i in seq_len(nrow(variant_manifest))) {
  row <- variant_manifest[i]
  path <- file.path(orientation_root, row$relative_path)
  require_hash(path, row$sha256)
  values <- fread(path)
  if (nrow(values) == 0L || !all(c(
    "snp", "position", "shared_posterior", "risk_to_expression_sign", "in_shared_95"
  ) %in% names(values))) stop("Malformed variant posterior: ", path)
  values[, orientation_uid := row$orientation_uid]
  variant_parts[[i]] <- values[as_flag(in_shared_95)]
}
variant <- rbindlist(variant_parts, use.names = TRUE, fill = TRUE)
variant <- merge(
  variant,
  worklist[, .(orientation_uid, chromosome)],
  by = "orientation_uid",
  all.x = TRUE,
  sort = FALSE
)
if (anyNA(variant$chromosome) || any(!is.finite(variant$shared_posterior))) {
  stop("Shared-95 variant assembly failed")
}

# Authoritative 1:1 hg19 -> hg38 conversion. A failed or multimapped position is
# retained as a failed routing row; it cannot overlap a peak.
chain_path <- file.path(project_root, "data/broadaway_eqtl/hg19ToHg38.over.chain")
if (!file.exists(chain_path)) stop("Missing hg19->hg38 chain")
chain <- import.chain(chain_path)
gr_hg19 <- GRanges(
  seqnames = paste0("chr", variant$chromosome),
  ranges = IRanges(start = as.integer(variant$position), width = 1L),
  strand = "*"
)
lifted <- liftOver(gr_hg19, chain)
variant[, liftover_n := lengths(lifted)]
variant[, `:=`(chr_hg38 = NA_character_, pos_hg38 = NA_integer_)]
unique_lift <- which(lengths(lifted) == 1L)
if (length(unique_lift) > 0L) {
  lifted_unique <- unlist(lifted[unique_lift])
  variant[unique_lift, `:=`(
    chr_hg38 = as.character(seqnames(lifted_unique)),
    pos_hg38 = start(lifted_unique)
  )]
}

load_peak_gr <- function(path) {
  peak <- unique(fread(path, select = 1:3))
  setnames(peak, c("chr", "start", "end"))
  GRanges(
    seqnames = peak$chr,
    ranges = IRanges(start = as.integer(peak$start) + 1L, end = as.integer(peak$end))
  )
}

# Revalidate every source path against the sealed reference before intersecting.
peak_gr <- list()
for (i in seq_len(nrow(access_manifest))) {
  source <- access_manifest[i]
  path <- file.path(project_root, source$source_path)
  require_hash(path, source$sha256)
  peak_gr[[source$source_id]] <- load_peak_gr(path)
}

routing <- merge(
  orientation,
  lineage[, .(
    gene_symbol,
    dominant_celltype,
    dominant_mean_cpm,
    dominant_detection_fraction_cpm1,
    dominant_to_second_ratio,
    lineage_expression_gate,
    lineage_specificity_gate,
    lineage_interpretation = interpretation
  )],
  by = "gene_symbol",
  all.x = TRUE,
  sort = FALSE
)
routing <- merge(
  routing,
  access_map[, .(
    expression_lineage,
    accessibility_lineage,
    mapping_status,
    eligible_for_source_lineage_gate,
    hepatocyte_independent_replication_required,
    mapping_note
  )],
  by.x = "dominant_celltype",
  by.y = "expression_lineage",
  all.x = TRUE,
  sort = FALSE
)

variant <- merge(
  variant,
  routing[, .(
    orientation_uid,
    gene_symbol,
    dominant_celltype,
    accessibility_lineage,
    mapping_status,
    orientation_gate_pass,
    pair_family_orientation_concordant
  )],
  by = "orientation_uid",
  all.x = TRUE,
  sort = FALSE
)
variant[, `:=`(
  primary_lineage_peak_overlap = FALSE,
  independent_hepatocyte_peak_overlap = FALSE
)]
valid_variant <- which(
  !is.na(variant$chr_hg38) &
    variant$mapping_status == "exact" &
    as_flag(variant$orientation_gate_pass) &
    as_flag(variant$pair_family_orientation_concordant)
)
for (index in valid_variant) {
  lineage_name <- variant$accessibility_lineage[[index]]
  primary_id <- paste0("GSE244832_", lineage_name)
  peaks <- peak_gr[[primary_id]]
  if (is.null(peaks)) next
  query <- GRanges(
    seqnames = variant$chr_hg38[[index]],
    ranges = IRanges(start = variant$pos_hg38[[index]], width = 1L)
  )
  variant$primary_lineage_peak_overlap[[index]] <- overlapsAny(query, peaks)
  if (identical(lineage_name, "Hepatocytes")) {
    variant$independent_hepatocyte_peak_overlap[[index]] <- overlapsAny(
      query,
      peak_gr[["GSE281367_Hepatocytes_IDR"]]
    )
  }
}
variant[, accessibility_variant_pass :=
  primary_lineage_peak_overlap &
  shared_posterior >= 0.05 &
  fifelse(
    accessibility_lineage == "Hepatocytes",
    independent_hepatocyte_peak_overlap,
    TRUE
  )
]

summary <- variant[, .(
  n_shared_95 = .N,
  n_liftover_1to1 = sum(liftover_n == 1L),
  n_primary_lineage_overlap = sum(primary_lineage_peak_overlap),
  primary_lineage_accessible_posterior = sum(
    shared_posterior[primary_lineage_peak_overlap]
  ),
  n_replicated_hepatocyte_overlap = sum(
    primary_lineage_peak_overlap & independent_hepatocyte_peak_overlap
  ),
  max_accessible_shared_posterior = suppressWarnings(max(
    shared_posterior[accessibility_variant_pass],
    na.rm = TRUE
  )),
  n_accessibility_variant_pass = sum(accessibility_variant_pass)
), by = orientation_uid]
summary[!is.finite(max_accessible_shared_posterior), max_accessible_shared_posterior := NA_real_]
routing <- merge(routing, summary, by = "orientation_uid", all.x = TRUE, sort = FALSE)
for (column in c(
  "n_shared_95", "n_liftover_1to1", "n_primary_lineage_overlap",
  "n_replicated_hepatocyte_overlap", "n_accessibility_variant_pass"
)) set(routing, which(is.na(routing[[column]])), column, 0L)
for (column in c("primary_lineage_accessible_posterior")) {
  set(routing, which(is.na(routing[[column]])), column, 0)
}
routing[, accessibility_evidence_class := fifelse(
  dominant_celltype == "Hepatocytes" & n_accessibility_variant_pass > 0L,
  "two_cohort_hepatocyte_accessibility",
  fifelse(
    n_accessibility_variant_pass > 0L,
    "single_cohort_exact_lineage_accessibility",
    "accessibility_gate_failed_or_unavailable"
  )
)]
routing[, routing_gate_pass :=
  as_flag(orientation_gate_pass) &
  as_flag(pair_family_orientation_concordant) &
  as_flag(lineage_expression_gate) &
  as_flag(lineage_specificity_gate) &
  mapping_status == "exact" &
  n_accessibility_variant_pass > 0L
]
routing[, routing_gate_reason := fifelse(
  routing_gate_pass,
  "pass_pending_guideability_and_exact_edit_feasibility",
  fifelse(
    !as_flag(orientation_gate_pass), "failed_pair_orientation",
    fifelse(
      !as_flag(pair_family_orientation_concordant), "failed_cross_gwas_orientation_concordance",
      fifelse(
        !as_flag(lineage_expression_gate), "failed_lineage_expression",
        fifelse(
          !as_flag(lineage_specificity_gate), "failed_unique_source_lineage",
          fifelse(
            mapping_status != "exact" | is.na(mapping_status), "failed_exact_accessibility_mapping",
            "failed_shared_signal_accessibility"
          )
        )
      )
    )
  )
)]
routing[, target_freeze_status := "prohibited_until_guideability_and_exact_edit_gate"]
variant[, target_freeze_status := "prohibited_until_guideability_and_exact_edit_gate"]

dir.create(candidate_root, recursive = TRUE)
variant_path <- file.path(candidate_root, "shared_signal_accessibility.tsv")
routing_path <- file.path(candidate_root, "lineage_accessibility_routing.tsv")
manifest_path <- file.path(candidate_root, "routing_input_manifest.tsv")
gate_path <- file.path(candidate_root, "routing_gate_status.tsv")
atomic_fwrite(variant, variant_path)
atomic_fwrite(routing, routing_path)

input_paths <- c(
  orientation_seal_path,
  lineage_seal_path,
  access_seal_path,
  chain_path
)
manifest <- data.table(
  role = c(
    "orientation_release", "lineage_reference", "accessibility_reference",
    "hg19_to_hg38_chain"
  ),
  source_path = sub(paste0("^", project_root, "/"), "", input_paths),
  size_bytes = file.info(input_paths)$size,
  sha256 = vapply(input_paths, sha256, character(1L))
)
atomic_fwrite(manifest, manifest_path)
gate <- data.table(
  gate = "lineage_accessibility_routing",
  n_prespecified_pairs = nrow(routing),
  n_pair_families = uniqueN(routing$pair_family_uid),
  n_physical_loci = uniqueN(routing$coarse_locus_uid),
  n_orientation_pair_pass = sum(as_flag(routing$orientation_gate_pass)),
  n_routing_pair_pass = sum(routing$routing_gate_pass),
  experimental_targets_frozen = FALSE,
  next_gate = "guideability_and_exact_edit_feasibility"
)
atomic_fwrite(gate, gate_path)

seal <- list(
  status = "lineage_accessibility_routing_complete_targets_not_selected",
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  candidate_id = candidate_id,
  n_prespecified_pairs = nrow(routing),
  n_pair_families = uniqueN(routing$pair_family_uid),
  n_physical_loci = uniqueN(routing$coarse_locus_uid),
  n_orientation_pair_gate_pass = sum(as_flag(routing$orientation_gate_pass)),
  n_routing_pair_gate_pass = sum(routing$routing_gate_pass),
  accessibility_variant_rule = paste(
    "shared-95 variant; shared posterior >=0.05; exact GSE244832 source-lineage peak;",
    "hepatocyte variants must also overlap the independent GSE281367 IDR peak set"
  ),
  source_lineage_rule = paste(
    "donor-collapsed dominant mean CPM>=1, donor detection>=0.25,",
    "dominant/second ratio>=2, and exact expression-to-ATAC lineage map"
  ),
  scientific_outcomes_inspected = FALSE,
  experimental_targets_frozen = FALSE,
  next_gate = "guideability and exact-edit feasibility; no ranking before that gate",
  output_sha256 = list(
    shared_signal_accessibility = sha256(variant_path),
    lineage_accessibility_routing = sha256(routing_path),
    routing_input_manifest = sha256(manifest_path),
    routing_gate_status = sha256(gate_path)
  )
)
seal_path <- file.path(candidate_root, "LINEAGE_ACCESSIBILITY_ROUTING_SEALED.json")
tmp <- tempfile(pattern = ".LINEAGE_ACCESSIBILITY_ROUTING_SEALED.", tmpdir = candidate_root)
writeLines(toJSON(seal, auto_unbox = TRUE, pretty = TRUE), tmp)
if (!file.rename(tmp, seal_path)) stop("Failed to seal routing release")
message(
  "[routing] completed ", nrow(routing), " gene-GWAS pairs; pass=",
  sum(routing$routing_gate_pass), "; targets not frozen"
)
