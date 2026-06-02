#!/usr/bin/env Rscript
# 11b_aggregate_mesusie.R
# Aggregate MESuSiE multi-ancestry fine-mapping results across all trait pairs
# into gene-level and variant-level summary tables.
#
# MESuSiE distinguishes shared (EUR_EAS), EUR-specific, and EAS-specific signals.
# pip_config columns give ancestry-conditional PIPs.
#
# Inputs:
#   results/mesusie/{trait_pair}/{locus_id}_mesusie.rds  — per-locus MESuSiE output
#   results/susiex/shared_loci.csv                       — locus metadata from Script 09
#   data/broadaway_eqtl/chr*_marginal_summary_results.tsv
#
# Outputs:
#   results/mesusie/mesusie_gene_summary.csv     — gene-level summary
#   results/mesusie/mesusie_variant_summary.csv  — variant-level PIPs
#
# Usage: Rscript 11b_aggregate_mesusie.R

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
})

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
FM_DIR <- Sys.getenv("FM_DIR",
  unset = file.path(Sys.getenv("MASLD_PROJECT_ROOT",
    unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"),
    "GWAS/finemapping"))
BASE_DIR <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

MESUSIE_RESULTS_SUFFIX <- Sys.getenv("MESUSIE_RESULTS_SUFFIX", unset = "")
MESUSIE_DIR     <- file.path(FM_DIR, paste0("results/mesusie", MESUSIE_RESULTS_SUFFIX))
SHARED_LOCI_CSV <- file.path(FM_DIR, "results/susiex/shared_loci.csv")
BROADAWAY_DIR   <- file.path(BASE_DIR, "data/broadaway_eqtl")
cat(sprintf("MESUSIE_DIR: %s (suffix='%s')\n", MESUSIE_DIR, MESUSIE_RESULTS_SUFFIX))

# ---------------------------------------------------------------------------
# 1. Read shared loci metadata
# ---------------------------------------------------------------------------
if (!file.exists(SHARED_LOCI_CSV)) {
  stop("shared_loci.csv not found — run 09_identify_shared_loci.R first.")
}
shared_loci <- read_csv(SHARED_LOCI_CSV, show_col_types = FALSE)
cat("Shared loci:", nrow(shared_loci), "\n")

# ---------------------------------------------------------------------------
# 2. Load all MESuSiE RDS files
# ---------------------------------------------------------------------------
rds_files <- list.files(MESUSIE_DIR, pattern = "_mesusie\\.rds$",
                         recursive = TRUE, full.names = TRUE)
cat("MESuSiE RDS files found:", length(rds_files), "\n")

if (length(rds_files) == 0) {
  message("No MESuSiE results found. Run 10b_run_mesusie.sh first.")
  quit(status = 0)
}

# ---------------------------------------------------------------------------
# 3. Parse each RDS into flat data.tables
# ---------------------------------------------------------------------------
variant_rows <- list()
locus_rows   <- list()

for (rds in rds_files) {
  obj <- tryCatch(readRDS(rds), error = function(e) NULL)
  if (is.null(obj)) next

  locus_id   <- obj$locus_id
  trait_pair <- obj$trait_pair
  snp_ids    <- obj$snp_ids
  pip        <- obj$pip           # overall PIP (marginal over all components)
  pip_config <- obj$pip_config    # matrix: nVariant x 3 (EUR, EAS, EUR_EAS)

  # Locus-level row
  n_shared  <- sum(obj$cs_category == "EUR_EAS", na.rm = TRUE)
  n_eur     <- sum(obj$cs_category == "EUR",     na.rm = TRUE)
  n_eas     <- sum(obj$cs_category == "EAS",     na.rm = TRUE)

  # Minimum CS size for each category
  cs_size_all <- vapply(obj$credible_sets, length, integer(1L))
  cs_size_shared <- cs_size_all[obj$cs_category == "EUR_EAS"]
  cs_size_eur    <- cs_size_all[obj$cs_category == "EUR"]
  cs_size_eas    <- cs_size_all[obj$cs_category == "EAS"]

  locus_rows[[length(locus_rows) + 1]] <- data.table(
    locus_id       = locus_id,
    trait_pair     = trait_pair,
    chr            = obj$chr,
    window_start   = obj$window[1],
    window_end     = obj$window[2],
    eur_gwas       = obj$eur_gwas,
    eas_gwas       = obj$eas_gwas,
    n_variants     = obj$n_variants,
    n_cs           = obj$n_cs,
    n_cs_shared    = n_shared,
    n_cs_eur       = n_eur,
    n_cs_eas       = n_eas,
    max_pip        = obj$max_pip,
    min_cs_size_shared = if (length(cs_size_shared) > 0) min(cs_size_shared) else NA_integer_,
    min_cs_size_eur    = if (length(cs_size_eur)    > 0) min(cs_size_eur)    else NA_integer_,
    min_cs_size_eas    = if (length(cs_size_eas)    > 0) min(cs_size_eas)    else NA_integer_,
    converged      = obj$converged
  )

  # Per-variant rows
  if (length(snp_ids) == 0 || length(pip) == 0) next

  vdt <- data.table(
    locus_id   = locus_id,
    trait_pair = trait_pair,
    variant_id = snp_ids,
    pip        = as.numeric(pip)
  )

  # Ancestry-conditional PIPs from pip_config
  if (!is.null(pip_config) && is.matrix(pip_config) &&
      nrow(pip_config) == length(snp_ids)) {
    if ("EUR"     %in% colnames(pip_config)) vdt$pip_eur     <- pip_config[, "EUR"]
    if ("EAS"     %in% colnames(pip_config)) vdt$pip_eas     <- pip_config[, "EAS"]
    if ("EUR_EAS" %in% colnames(pip_config)) vdt$pip_shared  <- pip_config[, "EUR_EAS"]
  }

  # CS membership: which CS does each variant belong to (index, or 0 = none)
  cs_member <- integer(length(snp_ids))
  cs_cat_vec <- character(length(snp_ids))
  if (!is.null(obj$credible_sets) && length(obj$credible_sets) > 0) {
    for (k in seq_along(obj$credible_sets)) {
      idx <- obj$credible_sets[[k]]
      cs_member[idx] <- k
      cs_cat_vec[idx] <- if (k <= length(obj$cs_category)) obj$cs_category[k] else ""
    }
  }
  vdt$cs_index    <- cs_member
  vdt$cs_category <- cs_cat_vec
  vdt$in_cs       <- cs_member > 0L

  variant_rows[[length(variant_rows) + 1]] <- vdt
}

all_variants <- rbindlist(variant_rows, fill = TRUE)
all_loci     <- rbindlist(locus_rows,   fill = TRUE)

cat("Total variants:", nrow(all_variants), "\n")
cat("Total loci:    ", nrow(all_loci),     "\n")
cat("Variants in CS:", sum(all_variants$in_cs, na.rm = TRUE), "\n")

# ---------------------------------------------------------------------------
# 4. Parse chr:pos from variant_id and load Broadaway eQTL gene map
# ---------------------------------------------------------------------------
cat("Loading Broadaway gene-variant map...\n")

parts <- strsplit(all_variants$variant_id, ":", fixed = TRUE)
all_variants[, `:=`(
  chr = as.integer(vapply(parts, `[[`, character(1), 1)),
  pos = as.integer(vapply(parts, `[[`, character(1), 2)),
  a1  = vapply(parts, `[[`, character(1), 3),
  a2  = vapply(parts, `[[`, character(1), 4)
)]

eqtl_files <- list.files(BROADAWAY_DIR,
  pattern = "_marginal_summary_results\\.tsv$", full.names = TRUE)
if (length(eqtl_files) == 0) stop("No Broadaway eQTL files found in ", BROADAWAY_DIR)

eqtl_map <- rbindlist(
  lapply(eqtl_files, function(f) {
    fread(f, select = c("GeneSymbol", "ENSG", "CHR", "POS", "NEA", "EA"))
  }),
  fill = TRUE
)
eqtl_map[, variant_id1 := paste(CHR, POS, EA,  NEA, sep = ":")]
eqtl_map[, variant_id2 := paste(CHR, POS, NEA, EA,  sep = ":")]
cat("  Broadaway gene-variant pairs:", nrow(eqtl_map), "\n")

# ---------------------------------------------------------------------------
# 5. Join variants to genes (allow allele flip)
# ---------------------------------------------------------------------------
fwd <- merge(all_variants, eqtl_map[, .(GeneSymbol, ENSG, variant_id1)],
             by.x = "variant_id", by.y = "variant_id1", all.x = FALSE)
rev <- merge(all_variants[!variant_id %in% fwd$variant_id],
             eqtl_map[, .(GeneSymbol, ENSG, variant_id2)],
             by.x = "variant_id", by.y = "variant_id2", all.x = FALSE)
all_joined <- rbindlist(list(fwd, rev), fill = TRUE)
cat("Variants with gene annotation:", nrow(all_joined), "\n")

# ---------------------------------------------------------------------------
# 6. Gene-level aggregation
# ---------------------------------------------------------------------------
gene_summary <- all_joined[, .(
  mesusie_max_pip        = max(pip, na.rm = TRUE),
  mesusie_max_pip_shared = max(pip_shared, na.rm = TRUE),
  mesusie_max_pip_eur    = max(pip_eur,    na.rm = TRUE),
  mesusie_max_pip_eas    = max(pip_eas,    na.rm = TRUE),
  mesusie_n_loci         = uniqueN(locus_id),
  mesusie_n_variants     = .N,
  mesusie_in_shared_cs   = any(in_cs & cs_category == "EUR_EAS", na.rm = TRUE),
  mesusie_in_eur_cs      = any(in_cs & cs_category == "EUR",     na.rm = TRUE),
  mesusie_in_eas_cs      = any(in_cs & cs_category == "EAS",     na.rm = TRUE),
  trait_pairs            = paste(sort(unique(trait_pair)), collapse = ";"),
  locus_ids              = paste(sort(unique(locus_id)),   collapse = ";")
), by = .(GeneSymbol, ENSG)]

# Replace -Inf from max(, na.rm=TRUE) on all-NA columns
inf_to_na <- function(x) { x[is.infinite(x)] <- NA_real_; x }
gene_summary[, mesusie_max_pip        := inf_to_na(mesusie_max_pip)]
gene_summary[, mesusie_max_pip_shared := inf_to_na(mesusie_max_pip_shared)]
gene_summary[, mesusie_max_pip_eur    := inf_to_na(mesusie_max_pip_eur)]
gene_summary[, mesusie_max_pip_eas    := inf_to_na(mesusie_max_pip_eas)]

# Add locus-level CS size info for the best locus per gene
best_locus_per_gene <- all_joined[in_cs == TRUE, .(
  locus_id  = locus_id[which.max(pip)],
  best_pip  = max(pip, na.rm = TRUE)
), by = .(GeneSymbol, ENSG)]

gene_cs <- merge(best_locus_per_gene[, .(GeneSymbol, ENSG, locus_id)],
                 all_loci[, .(locus_id, n_cs_shared, n_cs_eur, n_cs_eas,
                               min_cs_size_shared, min_cs_size_eur, min_cs_size_eas)],
                 by = "locus_id", all.x = TRUE)
gene_cs[, locus_id := NULL]

gene_summary <- merge(gene_summary, gene_cs, by = c("GeneSymbol", "ENSG"), all.x = TRUE)

setorder(gene_summary, -mesusie_max_pip)

cat("\n=== Top genes by MESuSiE max PIP ===\n")
print(head(gene_summary[, .(GeneSymbol, mesusie_max_pip, mesusie_max_pip_shared,
                             mesusie_in_shared_cs, n_cs_shared, min_cs_size_shared,
                             trait_pairs)], 20))

cat("\n=== Locus summary ===\n")
cat("Loci with converged fit: ", sum(all_loci$converged, na.rm = TRUE), "/",
    nrow(all_loci), "\n")
cat("Loci with >=1 shared CS: ", sum(all_loci$n_cs_shared > 0, na.rm = TRUE), "\n")
cat("Loci with >=1 EUR CS:    ", sum(all_loci$n_cs_eur   > 0, na.rm = TRUE), "\n")
cat("Loci with >=1 EAS CS:    ", sum(all_loci$n_cs_eas   > 0, na.rm = TRUE), "\n")
cat("Genes in shared CS:      ", sum(gene_summary$mesusie_in_shared_cs, na.rm = TRUE), "\n")
cat("Genes max PIP > 0.5:     ", sum(gene_summary$mesusie_max_pip > 0.5, na.rm = TRUE), "\n")
cat("Genes max PIP > 0.8:     ", sum(gene_summary$mesusie_max_pip > 0.8, na.rm = TRUE), "\n")

# ---------------------------------------------------------------------------
# 7. Write outputs
# ---------------------------------------------------------------------------
dir.create(MESUSIE_DIR, recursive = TRUE, showWarnings = FALSE)

out_gene    <- file.path(MESUSIE_DIR, "mesusie_gene_summary.csv")
out_variant <- file.path(MESUSIE_DIR, "mesusie_variant_summary.csv")
out_loci    <- file.path(MESUSIE_DIR, "mesusie_locus_summary.csv")

write_csv(as.data.frame(gene_summary),  out_gene)
write_csv(as.data.frame(all_variants),  out_variant)
write_csv(as.data.frame(all_loci),      out_loci)

cat("\nOutputs written:\n")
cat("  Gene-level:    ", out_gene,    "\n")
cat("  Variant-level: ", out_variant, "\n")
cat("  Locus-level:   ", out_loci,    "\n")
cat("  Genes with MESuSiE evidence:", nrow(gene_summary), "\n")
