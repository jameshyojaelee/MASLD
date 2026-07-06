#!/usr/bin/env Rscript
# 11b_aggregate_mesusie_mvp.R
# Aggregate within-MVP multi-ancestry MESuSiE results (EUR+EAS+AFR+AMR).
# N-ancestry generalization of 11b_aggregate_mesusie.R.
#
# Input: results/mesusie_mvp/{trait_pair}/{locus_id}_mesusie.rds
# Each RDS holds:
#   ancestries (chr vec), studies (named chr), pip (vec), pip_config (matrix
#   with one col per arm + N-1 shared mixtures), credible_sets, cs_category,
#   cs_purity, etc.
#
# Outputs:
#   results/mesusie_mvp/mesusie_gene_summary_mvp.csv
#   results/mesusie_mvp/mesusie_variant_summary_mvp.csv
#   results/mesusie_mvp/mesusie_locus_summary_mvp.csv

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(readr)
})

FM_DIR <- Sys.getenv("FM_DIR",
  unset = file.path(Sys.getenv("MASLD_PROJECT_ROOT",
    unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"),
    "GWAS/finemapping"))
BASE_DIR <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

MESUSIE_DIR     <- file.path(FM_DIR, "results/mesusie_mvp")
SHARED_LOCI_CSV <- file.path(FM_DIR, "results/susiex_mvp/shared_loci.csv")
BROADAWAY_DIR   <- file.path(BASE_DIR, "data/broadaway_eqtl")

if (!file.exists(SHARED_LOCI_CSV)) stop("shared_loci.csv missing")
shared_loci <- read_csv(SHARED_LOCI_CSV, show_col_types = FALSE)
cat("Shared loci:", nrow(shared_loci), "\n")

rds_files <- list.files(MESUSIE_DIR, pattern = "_mesusie\\.rds$",
                        recursive = TRUE, full.names = TRUE)
cat("MESuSiE 4-way RDS:", length(rds_files), "\n")
if (length(rds_files) == 0) { message("No RDS — exit"); quit(status = 0) }

variant_rows <- list(); locus_rows <- list()

for (rds in rds_files) {
  obj <- tryCatch(readRDS(rds), error = function(e) NULL)
  if (is.null(obj)) next
  lid <- obj$locus_id; tp <- obj$trait_pair
  ancs <- obj$ancestries
  snp_ids <- obj$snp_ids; pip <- obj$pip; pip_config <- obj$pip_config

  # Locus row: capture per-ancestry n_cs and per-category min CS size
  cs_cat <- if (!is.null(obj$cs_category)) obj$cs_category else character(0)
  cs_size_all <- if (!is.null(obj$credible_sets)) vapply(obj$credible_sets, length, integer(1)) else integer(0)
  cat_table <- table(cs_cat)

  loc_row <- data.table(
    locus_id   = lid,
    trait_pair = tp,
    chr        = obj$chr,
    window_start = obj$window[1],
    window_end   = obj$window[2],
    n_arms     = length(ancs),
    arms       = paste(ancs, collapse = ","),
    studies    = paste(obj$studies, collapse = ","),
    n_variants = obj$n_variants,
    n_cs       = obj$n_cs,
    max_pip    = obj$max_pip,
    converged  = obj$converged
  )
  # Add count + min CS size per category as dynamic columns
  for (cat_name in names(cat_table)) {
    loc_row[[paste0("n_cs_", cat_name)]] <- as.integer(cat_table[[cat_name]])
    loc_row[[paste0("min_cs_size_", cat_name)]] <- min(cs_size_all[cs_cat == cat_name])
  }
  locus_rows[[length(locus_rows) + 1]] <- loc_row

  if (length(snp_ids) == 0 || length(pip) == 0) next

  vdt <- data.table(locus_id = lid, trait_pair = tp,
                    variant_id = snp_ids, pip = as.numeric(pip))

  # Per-component PIP from pip_config matrix (one column per ancestry + shared mixtures)
  if (!is.null(pip_config) && is.matrix(pip_config) && nrow(pip_config) == length(snp_ids)) {
    for (col in colnames(pip_config)) {
      vdt[[paste0("pip_", col)]] <- pip_config[, col]
    }
  }

  # CS membership
  cs_member <- integer(length(snp_ids))
  cs_cat_vec <- character(length(snp_ids))
  if (!is.null(obj$credible_sets) && length(obj$credible_sets) > 0) {
    for (k in seq_along(obj$credible_sets)) {
      idx <- obj$credible_sets[[k]]
      cs_member[idx] <- k
      cs_cat_vec[idx] <- if (k <= length(cs_cat)) cs_cat[k] else ""
    }
  }
  vdt$cs_index    <- cs_member
  vdt$cs_category <- cs_cat_vec
  vdt$in_cs       <- cs_member > 0L
  variant_rows[[length(variant_rows) + 1]] <- vdt
}

all_variants <- rbindlist(variant_rows, fill = TRUE)
all_loci     <- rbindlist(locus_rows,   fill = TRUE)
cat("Variants total:", nrow(all_variants), "\n")
cat("Variants in CS:", sum(all_variants$in_cs, na.rm = TRUE), "\n")
cat("Loci:          ", nrow(all_loci), "\n")

# Parse chr/pos from SNP id
parts <- strsplit(all_variants$variant_id, ":", fixed = TRUE)
all_variants[, `:=`(chr = as.integer(vapply(parts, `[[`, character(1), 1)),
                    pos = as.integer(vapply(parts, `[[`, character(1), 2)),
                    a1  = vapply(parts, `[[`, character(1), 3),
                    a2  = vapply(parts, `[[`, character(1), 4))]

# Broadaway gene-variant map
cat("Loading Broadaway eQTL...\n")
eqtl_files <- list.files(BROADAWAY_DIR, pattern = "_marginal_summary_results\\.tsv$", full.names = TRUE)
eqtl_map <- rbindlist(lapply(eqtl_files, function(f)
  fread(f, select = c("GeneSymbol", "ENSG", "CHR", "POS", "NEA", "EA"))), fill = TRUE)
eqtl_map[, variant_id1 := paste(CHR, POS, EA, NEA, sep = ":")]
eqtl_map[, variant_id2 := paste(CHR, POS, NEA, EA, sep = ":")]

fwd <- merge(all_variants, eqtl_map[, .(GeneSymbol, ENSG, variant_id1)],
             by.x = "variant_id", by.y = "variant_id1", all.x = FALSE)
rev <- merge(all_variants[!variant_id %in% fwd$variant_id],
             eqtl_map[, .(GeneSymbol, ENSG, variant_id2)],
             by.x = "variant_id", by.y = "variant_id2", all.x = FALSE)
all_joined <- rbindlist(list(fwd, rev), fill = TRUE)
cat("Variants with gene annotation:", nrow(all_joined), "\n")

# Gene-level summary
inf_to_na <- function(x) { x[is.infinite(x)] <- NA_real_; x }

# Discover all pip_<...> columns dynamically
pip_cols <- grep("^pip_", names(all_joined), value = TRUE)

gene_summary <- all_joined[, .(
  mesusie_max_pip    = max(pip, na.rm = TRUE),
  mesusie_n_loci     = uniqueN(locus_id),
  mesusie_n_variants = .N,
  mesusie_in_any_cs  = any(in_cs, na.rm = TRUE),
  cs_categories      = paste(sort(unique(cs_category[in_cs])), collapse = ";"),
  trait_pairs        = paste(sort(unique(trait_pair)), collapse = ";"),
  locus_ids          = paste(sort(unique(locus_id)),   collapse = ";")
), by = .(GeneSymbol, ENSG)]

# Add max per pip_<config> column
for (pc in pip_cols) {
  agg <- all_joined[, .(maxv = max(get(pc), na.rm = TRUE)), by = .(GeneSymbol, ENSG)]
  agg[, maxv := inf_to_na(maxv)]
  setnames(agg, "maxv", paste0("max_", pc))
  gene_summary <- merge(gene_summary, agg, by = c("GeneSymbol", "ENSG"), all.x = TRUE)
}
gene_summary[, mesusie_max_pip := inf_to_na(mesusie_max_pip)]
setorder(gene_summary, -mesusie_max_pip)

cat("\n=== Top genes by MESuSiE 4-way max PIP ===\n")
print(head(gene_summary[, .(GeneSymbol, mesusie_max_pip, mesusie_in_any_cs,
                             cs_categories, trait_pairs)], 20))

cat(sprintf("\nLoci converged: %d/%d  Genes PIP>0.5: %d  Genes PIP>0.8: %d\n",
            sum(all_loci$converged, na.rm = TRUE), nrow(all_loci),
            sum(gene_summary$mesusie_max_pip > 0.5, na.rm = TRUE),
            sum(gene_summary$mesusie_max_pip > 0.8, na.rm = TRUE)))

dir.create(MESUSIE_DIR, recursive = TRUE, showWarnings = FALSE)
write_csv(as.data.frame(gene_summary), file.path(MESUSIE_DIR, "mesusie_gene_summary_mvp.csv"))
write_csv(as.data.frame(all_variants), file.path(MESUSIE_DIR, "mesusie_variant_summary_mvp.csv"))
write_csv(as.data.frame(all_loci),     file.path(MESUSIE_DIR, "mesusie_locus_summary_mvp.csv"))
cat("Outputs in:", MESUSIE_DIR, "\n")
