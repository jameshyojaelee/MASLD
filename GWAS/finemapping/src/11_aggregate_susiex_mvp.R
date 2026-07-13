#!/usr/bin/env Rscript
# 11_aggregate_susiex_mvp.R
# Aggregate within-MVP multi-ancestry SuSiEX results (EUR+EAS+AFR+AMR, 1KG LD).
# Outputs per-variant + per-gene summary CSVs analogous to the 2-way aggregators.
#
# Inputs:
#   results/susiex_mvp/{trait_pair}/{locus_id}.{snp,cs,summary}
#   results/susiex_mvp/shared_loci.csv                  — locus metadata (LD/ancestry agnostic)
#   data/broadaway_eqtl/chr*_marginal_summary_results.tsv — gene-variant map
#
# Outputs:
#   results/susiex_mvp/susiex_gene_summary_mvp.csv     — per-gene max PIP, CS counts
#   results/susiex_mvp/susiex_variant_summary_mvp.csv  — per-variant PIPs
#   results/susiex_mvp/susiex_cs_mvp.csv               — credible sets w/ per-pop posteriors
#
# Usage: Rscript 11_aggregate_susiex_mvp.R

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(readr); library(purrr); library(tidyr)
})

FM_DIR   <- Sys.getenv("FM_DIR",
  unset = file.path(Sys.getenv("MASLD_PROJECT_ROOT",
    unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"),
    "GWAS/finemapping"))
BASE_DIR <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

SUSIEX_DIR      <- file.path(FM_DIR, "results/susiex_mvp")
SHARED_LOCI_CSV <- file.path(FM_DIR, "results/susiex_mvp/shared_loci.csv")
BROADAWAY_DIR   <- file.path(BASE_DIR, "data/broadaway_eqtl")

if (!file.exists(SHARED_LOCI_CSV)) stop("shared_loci.csv missing: ", SHARED_LOCI_CSV)
shared_loci <- read_csv(SHARED_LOCI_CSV, show_col_types = FALSE)
cat("Shared loci:", nrow(shared_loci), "\n")

# --- Load .snp (per-variant PIPs)
load_snp <- function(snp_file, lid, tp) {
  if (!file.exists(snp_file)) return(NULL)
  dt <- tryCatch(fread(snp_file), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  pip_cols <- grep("^PIP\\(CS", names(dt), value = TRUE)
  ovrl_col <- grep("^OVRL_PIP$", names(dt), value = TRUE)
  if (length(ovrl_col) == 0 && length(pip_cols) > 0) {
    dt[, OVRL_PIP := apply(.SD, 1, max), .SDcols = pip_cols]
    ovrl_col <- "OVRL_PIP"
  }
  if (length(ovrl_col) == 0) return(NULL)
  out <- dt[, .(BP, SNP, OVRL_PIP = get(ovrl_col[1]))]
  out[, `:=`(locus_id = lid, trait_pair = tp,
             chr = as.integer(sub(":.*", "", SNP)))]
  out
}

# --- Load .cs (one row per variant in a CS; CS_ID groups variants)
load_cs <- function(cs_file, lid, tp) {
  if (!file.exists(cs_file)) return(NULL)
  dt <- tryCatch(fread(cs_file), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  dt[, `:=`(locus_id = lid, trait_pair = tp)]
  dt
}

trait_dirs <- list.dirs(SUSIEX_DIR, recursive = FALSE, full.names = TRUE)
snp_list <- list(); cs_list <- list()
for (tdir in trait_dirs) {
  tp <- basename(tdir)
  for (f in list.files(tdir, pattern = "\\.snp$", full.names = TRUE)) {
    lid <- tools::file_path_sans_ext(basename(f))
    d <- load_snp(f, lid, tp); if (!is.null(d)) snp_list[[length(snp_list)+1]] <- d
  }
  for (f in list.files(tdir, pattern = "\\.cs$", full.names = TRUE)) {
    lid <- tools::file_path_sans_ext(basename(f))
    d <- load_cs(f, lid, tp); if (!is.null(d)) cs_list[[length(cs_list)+1]] <- d
  }
}
if (length(snp_list) == 0) { message("No .snp files in ", SUSIEX_DIR); quit(status = 0) }

all_snp <- rbindlist(snp_list, fill = TRUE)
all_cs  <- if (length(cs_list) > 0) rbindlist(cs_list, fill = TRUE) else NULL
cat("Variants with PIP:", nrow(all_snp), "\n")
cat("Variants in CS:   ", if (!is.null(all_cs)) nrow(all_cs) else 0, "\n")

# Parse SNP id → chr/pos/a1/a2
parts <- strsplit(all_snp$SNP, ":", fixed = TRUE)
all_snp[, `:=`(chr = as.integer(sapply(parts, `[`, 1)),
               pos = as.integer(sapply(parts, `[`, 2)),
               a1  = sapply(parts, `[`, 3),
               a2  = sapply(parts, `[`, 4))]
all_snp[, variant_id := paste(chr, pos, a1, a2, sep = ":")]

# --- Broadaway gene-variant map (allele-flip aware)
cat("Loading Broadaway eQTL ...\n")
eqtl_files <- list.files(BROADAWAY_DIR, pattern = "_marginal_summary_results\\.tsv$", full.names = TRUE)
if (length(eqtl_files) == 0) stop("No Broadaway files in ", BROADAWAY_DIR)
eqtl_map <- rbindlist(lapply(eqtl_files, function(f)
  fread(f, select = c("GeneSymbol", "ENSG", "CHR", "POS", "NEA", "EA"))), fill = TRUE)
eqtl_map[, variant_id1 := paste(CHR, POS, EA, NEA, sep = ":")]
eqtl_map[, variant_id2 := paste(CHR, POS, NEA, EA, sep = ":")]

joined <- merge(all_snp, eqtl_map[, .(GeneSymbol, ENSG, variant_id1)],
                by.x = "variant_id", by.y = "variant_id1", all.x = FALSE)
joined_flip <- merge(all_snp[!variant_id %in% joined$variant_id],
                     eqtl_map[, .(GeneSymbol, ENSG, variant_id2)],
                     by.x = "variant_id", by.y = "variant_id2", all.x = FALSE)
all_joined <- rbindlist(list(joined, joined_flip), fill = TRUE)
cat("Variants with gene annotation:", nrow(all_joined), "\n")

# --- CS sizes per locus
if (!is.null(all_cs)) {
  cs_sizes <- all_cs[, .(susiex_cs_size_joint = uniqueN(SNP),
                         susiex_n_cs = uniqueN(CS_ID)), by = locus_id]
} else {
  cs_sizes <- data.table(locus_id = character(), susiex_cs_size_joint = integer(),
                          susiex_n_cs = integer())
}

# --- Gene-level summary
gene_summary <- all_joined[, .(
  susiex_max_pip    = max(OVRL_PIP, na.rm = TRUE),
  susiex_n_loci     = uniqueN(locus_id),
  susiex_n_variants = .N,
  trait_pairs       = paste(sort(unique(trait_pair)), collapse = ";"),
  locus_ids         = paste(sort(unique(locus_id)), collapse = ";")
), by = .(GeneSymbol, ENSG)]

gene_cs <- merge(all_joined[OVRL_PIP > 0.1], cs_sizes, by = "locus_id", all.x = TRUE)
gene_cs_agg <- gene_cs[, .(susiex_n_cs = max(susiex_n_cs, na.rm = TRUE),
                            susiex_cs_size_joint = min(susiex_cs_size_joint, na.rm = TRUE)),
                        by = .(GeneSymbol, ENSG)]
gene_cs_agg[is.infinite(susiex_n_cs), susiex_n_cs := NA_integer_]
gene_cs_agg[is.infinite(susiex_cs_size_joint), susiex_cs_size_joint := NA_integer_]
gene_summary <- merge(gene_summary, gene_cs_agg, by = c("GeneSymbol", "ENSG"), all.x = TRUE)
setorder(gene_summary, -susiex_max_pip)

cat("\n=== Top genes by 4-way SuSiEX PIP ===\n")
print(head(gene_summary[, .(GeneSymbol, susiex_max_pip, susiex_n_cs,
                             susiex_cs_size_joint, trait_pairs)], 20))

# --- Write outputs
out_gene    <- file.path(SUSIEX_DIR, "susiex_gene_summary_mvp.csv")
out_variant <- file.path(SUSIEX_DIR, "susiex_variant_summary_mvp.csv")
out_cs      <- file.path(SUSIEX_DIR, "susiex_cs_mvp.csv")
write_csv(as.data.frame(gene_summary), out_gene)
write_csv(as.data.frame(all_snp),      out_variant)
if (!is.null(all_cs)) write_csv(as.data.frame(all_cs), out_cs)
cat("\nOutputs:\n  ", out_gene, "\n  ", out_variant, "\n  ", out_cs, "\n")
cat("Genes with PIP > 0.5:", nrow(gene_summary[susiex_max_pip > 0.5]), "\n")
cat("Loci aggregated:    ", uniqueN(all_snp$locus_id), "\n")
