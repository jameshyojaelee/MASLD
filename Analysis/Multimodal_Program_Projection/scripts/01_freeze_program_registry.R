#!/usr/bin/env Rscript

# Freeze the exact donor-level Hotspot program contract used by Figure 3.
# This script is intentionally assay-blind: no ATAC, spatial, or proteomic data
# are read while program identities, membership, weights, and display order are
# fixed.

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})

BASE <- Sys.getenv(
  "MASLD_PROJECT_ROOT",
  "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
ROOT <- file.path(BASE, "Analysis/Multimodal_Program_Projection")
OUT <- file.path(ROOT, "results")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

sha256 <- function(path) digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
CONTRACT_FILE <- file.path(ROOT, "config/fig4_input_contract.tsv")
if (!file.exists(CONTRACT_FILE)) stop("Missing versioned Figure 4 input contract: ", CONTRACT_FILE)
contract <- fread(CONTRACT_FILE)
stopifnot(all(c("role", "relative_path", "expected_sha256") %in% names(contract)))
fig3_contract <- contract[role == "fig3_input"]
fig3_contract[, path := file.path(BASE, relative_path)]
if (any(!file.exists(fig3_contract$path))) {
  stop(
    "Missing pinned Figure 3 input(s): ",
    paste(fig3_contract[!file.exists(path), relative_path], collapse = ", ")
  )
}
fig3_contract[, observed_sha256 := vapply(path, sha256, character(1))]
fig3_contract[, contract_match := observed_sha256 == expected_sha256]
if (any(!fig3_contract$contract_match)) {
  bad <- fig3_contract[contract_match == FALSE]
  stop(
    "Figure 3 input drifted from the versioned Figure 4 contract: ",
    paste(bad$relative_path, collapse = ", "),
    ". Review the upstream change and deliberately update config/fig4_input_contract.tsv before rebuilding."
  )
}

FIG3 <- file.path(
  BASE, "figures/main/fig3_RNAseq/panels/data/fig3g_module_matrix.csv"
)
DONOR <- file.path(
  BASE,
  "Analysis/SingleCell/results_gpu_v2/hotspot_modules/donor_collapse/all_modules_donor.tsv"
)
ATLAS <- file.path(
  BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"
)
MODULE_ROOT <- file.path(
  BASE, "Analysis/SingleCell/results_gpu_v2/hotspot_modules"
)

inputs <- c(fig3_matrix = FIG3, donor_modules = DONOR, atlas = ATLAS)
missing <- names(inputs)[!file.exists(inputs)]
if (length(missing)) stop("Missing freeze input(s): ", paste(missing, collapse = ", "))

fig3 <- fread(FIG3)
donor <- fread(DONOR)
stopifnot(
  nrow(fig3) == 22L,
  all(c("cell_type", "module", "name", "beta", "q", "dir") %in% names(fig3)),
  all(c("cell_type", "module", "disease_stage_q", "stability_fail") %in% names(donor))
)

fig3[, `:=`(
  module = as.integer(module),
  program_id = sprintf("%s::%d", cell_type, as.integer(module)),
  display_order = .I
)]
donor[, `:=`(
  module = as.integer(module),
  program_id = sprintf("%s::%d", cell_type, as.integer(module))
)]

disease22 <- donor[is.finite(disease_stage_q) & disease_stage_q < 0.05]
if (nrow(disease22) != 22L) {
  stop("Expected 22 donor-level disease modules at q<0.05; observed ", nrow(disease22))
}
if (!setequal(fig3$program_id, disease22$program_id)) {
  stop("Figure 3 rows do not equal the donor-level q<0.05 Hotspot set")
}

registry <- merge(
  fig3,
  donor,
  by = c("cell_type", "module", "program_id"),
  all.x = TRUE,
  sort = FALSE,
  suffixes = c("_fig3", "_donor")
)
setorder(registry, display_order)
registry[, lineage_order := match(
  cell_type,
  c("hepatocytes", "fibroblasts", "macrophages", "cholangiocytes")
)]
registry[, direction := fifelse(beta >= 0, "up", "down")]
if (any(registry$direction != registry$dir)) {
  stop("Figure 3 direction labels disagree with donor-level coefficients")
}

# Resolve Ensembl IDs without altering source membership. The atlas is preferred
# because it is the paper's canonical symbol map; unmapped non-Ensembl entries
# are already gene symbols.
atlas <- fread(ATLAS, select = c("ensembl_id", "human_symbol", "gene_biotype"))
atlas[, ensembl_base := sub("\\..*$", "", ensembl_id)]
atlas <- atlas[!is.na(ensembl_base) & ensembl_base != ""]
atlas[, symbol_missing := is.na(human_symbol) | human_symbol == ""]
setorder(atlas, ensembl_base, symbol_missing, human_symbol)
atlas <- atlas[!duplicated(ensembl_base)]
symbol_map <- setNames(atlas$human_symbol, atlas$ensembl_base)
biotype_map <- setNames(atlas$gene_biotype, atlas$ensembl_base)

members <- list()
for (ct in unique(registry$cell_type)) {
  f <- file.path(MODULE_ROOT, ct, "module_genes.tsv")
  if (!file.exists(f)) stop("Missing module membership file: ", f)
  inputs[paste0("membership_", ct)] <- f
  x <- fread(f)
  x[, module := as.integer(module)]
  x <- x[module %in% registry[cell_type == ct, module]]
  x[, `:=`(
    cell_type = ct,
    source_gene = as.character(gene),
    ensembl_base = fifelse(grepl("^ENSG", gene), sub("\\..*$", "", gene), NA_character_)
  )]
  x[, gene_symbol := fifelse(
    !is.na(ensembl_base), unname(symbol_map[ensembl_base]), source_gene
  )]
  x[, gene_biotype := fifelse(
    !is.na(ensembl_base), unname(biotype_map[ensembl_base]), NA_character_
  )]
  x[, mapped_symbol := !is.na(gene_symbol) & gene_symbol != ""]
  x[, program_id := sprintf("%s::%d", cell_type, module)]
  x[, weight := as.numeric(weight)]
  if (any(!is.finite(x$weight) | x$weight <= 0)) {
    stop("Non-positive or non-finite Hotspot weight in ", f)
  }
  x[, original_l1_weight := weight / sum(abs(weight)), by = program_id]
  members[[ct]] <- x[, .(
    program_id, cell_type, module, source_gene, ensembl_base, gene_symbol,
    gene_biotype, mapped_symbol, weight, original_l1_weight
  )]
}
membership <- rbindlist(members, use.names = TRUE)

member_summary <- membership[, .(
  n_source_genes = .N,
  n_symbol_mapped = sum(mapped_symbol),
  mapped_l1_weight = sum(original_l1_weight[mapped_symbol])
), by = program_id]
registry <- merge(registry, member_summary, by = "program_id", all.x = TRUE, sort = FALSE)
setorder(registry, display_order)
if (anyDuplicated(registry$program_id) || any(is.na(registry$n_source_genes))) {
  stop("Program registry is not one-to-one with membership")
}

registry_out <- registry[, .(
  program_id,
  display_order,
  lineage_order,
  cell_type,
  module,
  program_name = name,
  fig3_beta = beta,
  fig3_se = SE,
  fig3_q = q,
  fig3_direction = direction,
  stability_score,
  stability_fail,
  gse244832_loo_jaccard,
  n_source_genes,
  n_symbol_mapped,
  mapped_l1_weight
)]

registry_file <- file.path(OUT, "frozen_programs.tsv")
membership_file <- file.path(OUT, "frozen_program_membership.tsv")
fwrite(registry_out, registry_file, sep = "\t", quote = FALSE)
fwrite(membership, membership_file, sep = "\t", quote = FALSE)

manifest <- data.table(
  role = c("contract", rep("input", length(inputs)), "output", "output"),
  path = c(CONTRACT_FILE, unname(inputs), registry_file, membership_file)
)
manifest[, relative_path := sub(paste0("^", BASE, "/?"), "", path)]
manifest[, `:=`(
  bytes = file.info(path)$size,
  sha256 = vapply(path, sha256, character(1))
)]
manifest[, expected_sha256 := fig3_contract$expected_sha256[match(relative_path, fig3_contract$relative_path)]]
manifest[, contract_match := fifelse(
  is.na(expected_sha256),
  NA,
  sha256 == expected_sha256
)]
manifest[, frozen_at := format(Sys.time(), tz = "UTC", usetz = TRUE)]
fwrite(
  manifest[, .(role, relative_path, bytes, sha256, expected_sha256, contract_match, frozen_at)],
  file.path(OUT, "freeze_manifest.tsv"),
  sep = "\t",
  quote = FALSE
)

cat("[freeze] 22 donor-level programs frozen\n")
cat("[freeze] lineage counts:\n")
print(registry_out[, .N, by = cell_type])
cat("[freeze] stability_fail:", sum(registry_out$stability_fail), "/22\n")
cat("[freeze] wrote:", registry_file, "\n")
cat("[freeze] wrote:", membership_file, "\n")
