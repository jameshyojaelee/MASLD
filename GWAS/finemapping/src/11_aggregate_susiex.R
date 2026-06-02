#!/usr/bin/env Rscript
# 11_aggregate_susiex.R
# Aggregate SuSiEX results across all trait pairs into a gene-level summary,
# and compare joint credible set sizes against single-ancestry SuSiE results.
#
# Inputs:
#   results/susiex/{trait_pair}/{locus_id}.cs   — SuSiEX per-variant CS output
#   results/susiex/{trait_pair}/{locus_id}.snp  — SuSiEX per-variant PIP output
#   results/susiex/shared_loci.csv              — locus metadata from script 09
#   results/combined_finemapping.csv            — single-ancestry SuSiE/CARMA results
#   data/broadaway_eqtl/chr*_marginal_summary_results.tsv — gene-variant mappings
#
# Outputs:
#   results/susiex/susiex_gene_summary.csv     — gene-level max PIP, CS counts,
#                                                 CS size comparison
#   results/susiex/susiex_variant_summary.csv  — variant-level PIP table
#
# Usage: Rscript 11_aggregate_susiex.R
#   Override root: MASLD_PROJECT_ROOT=/path Rscript 11_aggregate_susiex.R

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
  library(purrr)
  library(tidyr)
})

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
FM_DIR   <- Sys.getenv("FM_DIR",
  unset = file.path(Sys.getenv("MASLD_PROJECT_ROOT",
    unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"),
    "GWAS/finemapping"))
BASE_DIR <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")

SUSIEX_DIR      <- file.path(FM_DIR, "results/susiex")
COMBINED_FM     <- file.path(FM_DIR, "results/combined_finemapping.csv")
BROADAWAY_DIR   <- file.path(BASE_DIR, "data/broadaway_eqtl")
SHARED_LOCI_CSV <- file.path(SUSIEX_DIR, "shared_loci.csv")

# ---------------------------------------------------------------------------
# 1. Read shared loci metadata
# ---------------------------------------------------------------------------
if (!file.exists(SHARED_LOCI_CSV)) {
  stop("shared_loci.csv not found. Run 09_identify_shared_loci.R first.")
}
shared_loci <- read_csv(SHARED_LOCI_CSV, show_col_types = FALSE)
cat("Shared loci:", nrow(shared_loci), "\n")
if (nrow(shared_loci) == 0) {
  message("No shared loci found — nothing to aggregate.")
  quit(status = 0)
}

# ---------------------------------------------------------------------------
# 2. Load SuSiEX per-variant output files
# ---------------------------------------------------------------------------
# .snp files contain per-variant PIP (all variants used in fine-mapping)
# .cs  files contain variants assigned to credible sets

load_susiex_snp <- function(snp_file, locus_id, trait_pair) {
  # Header: BP  SNP  PIP(CS1)  PIP(CS2) ...  OVRL_PIP
  # Also may contain: LogBF(CS1,Pop1) etc.
  # We extract: BP, SNP, and the OVRL_PIP column
  if (!file.exists(snp_file)) return(NULL)
  content <- readLines(snp_file, n = 2)
  if (length(content) == 0 || grepl("^NULL|^FAIL", content[1])) return(NULL)
  dt <- tryCatch(fread(snp_file), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  nms <- names(dt)
  # Identify PIP columns (named PIP(CS*))
  pip_cols <- grep("^PIP\\(CS", nms, value = TRUE)
  ovrl_col <- grep("^OVRL_PIP$", nms, value = TRUE)
  if (length(ovrl_col) == 0 && length(pip_cols) > 0) {
    # Compute OVRL_PIP as max across CS PIPs if not present
    dt[, OVRL_PIP := apply(.SD, 1, max), .SDcols = pip_cols]
    ovrl_col <- "OVRL_PIP"
  }
  if (length(ovrl_col) == 0) return(NULL)
  result <- dt[, .(BP, SNP, OVRL_PIP = get(ovrl_col[1]))]
  result[, `:=`(locus_id = locus_id, trait_pair = trait_pair,
                chr = as.integer(sub(":.*", "", SNP)))]
  result
}

load_susiex_cs <- function(cs_file, locus_id, trait_pair) {
  # Header: CS_ID  SNP  BP  REF_ALLELE  ALT_ALLELE  REF_FRQ  BETA  SE  -LOG10P  CS_PIP  OVRL_PIP
  if (!file.exists(cs_file)) return(NULL)
  content <- readLines(cs_file, n = 2)
  if (length(content) == 0 || grepl("^NULL|^FAIL", content[1])) return(NULL)
  dt <- tryCatch(fread(cs_file), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  dt[, `:=`(locus_id = locus_id, trait_pair = trait_pair)]
  dt
}

# Discover result directories
trait_dirs <- list.dirs(SUSIEX_DIR, recursive = FALSE, full.names = TRUE)
trait_dirs <- trait_dirs[!grepl("/$", trait_dirs)]  # remove trailing slash dirs

snp_list <- list()
cs_list  <- list()

for (tdir in trait_dirs) {
  tp <- basename(tdir)
  snp_files <- list.files(tdir, pattern = "\\.snp$", full.names = TRUE)
  cs_files  <- list.files(tdir, pattern = "\\.cs$",  full.names = TRUE)

  for (f in snp_files) {
    lid <- tools::file_path_sans_ext(basename(f))
    dat <- load_susiex_snp(f, lid, tp)
    if (!is.null(dat)) snp_list[[length(snp_list) + 1]] <- dat
  }
  for (f in cs_files) {
    lid <- tools::file_path_sans_ext(basename(f))
    dat <- load_susiex_cs(f, lid, tp)
    if (!is.null(dat)) cs_list[[length(cs_list) + 1]] <- dat
  }
}

if (length(snp_list) == 0) {
  message("No SuSiEX .snp files found. Run 10_run_susiex.sh first.")
  quit(status = 0)
}

all_snp <- rbindlist(snp_list, fill = TRUE)
all_cs  <- if (length(cs_list) > 0) rbindlist(cs_list, fill = TRUE) else NULL

cat("Variants with PIP:", nrow(all_snp), "\n")
cat("Variants in CS:",    if (!is.null(all_cs)) nrow(all_cs) else 0, "\n")

# ---------------------------------------------------------------------------
# 3. Variant-level summary table
# ---------------------------------------------------------------------------
# Parse chr and pos from SNP ID (format: CHR:POS:A1:A2)
parse_variant <- function(dt) {
  parts <- strsplit(dt$SNP, ":", fixed = TRUE)
  dt$chr  <- as.integer(sapply(parts, `[[`, 1))
  dt$pos  <- as.integer(sapply(parts, `[[`, 2))
  dt$a1   <- sapply(parts, `[[`, 3)
  dt$a2   <- sapply(parts, `[[`, 4)
  dt$variant_id <- paste(dt$chr, dt$pos, dt$a1, dt$a2, sep = ":")
  dt
}

all_snp <- parse_variant(as.data.frame(all_snp))

# ---------------------------------------------------------------------------
# 4. Load Broadaway eQTL gene-variant map
# ---------------------------------------------------------------------------
cat("Loading Broadaway eQTL gene annotations ...\n")
eqtl_files <- list.files(BROADAWAY_DIR, pattern = "_marginal_summary_results\\.tsv$",
                          full.names = TRUE)
if (length(eqtl_files) == 0) {
  stop("No Broadaway eQTL files found in ", BROADAWAY_DIR)
}

eqtl_map <- rbindlist(
  lapply(eqtl_files, function(f) {
    dt <- fread(f, select = c("GeneSymbol", "ENSG", "CHR", "POS", "NEA", "EA"))
    dt
  }),
  fill = TRUE
)
# Build variant_id to match all_snp: CHR:POS:EA:NEA or CHR:POS:NEA:EA
eqtl_map[, variant_id1 := paste(CHR, POS, EA, NEA, sep = ":")]
eqtl_map[, variant_id2 := paste(CHR, POS, NEA, EA, sep = ":")]
cat("  Broadaway gene-variant pairs:", nrow(eqtl_map), "\n")

# ---------------------------------------------------------------------------
# 5. Join SuSiEX variants to gene annotations (allow allele flip)
# ---------------------------------------------------------------------------
all_snp_dt <- as.data.table(all_snp)

joined <- merge(
  all_snp_dt,
  eqtl_map[, .(GeneSymbol, ENSG, variant_id1, variant_id2)],
  by.x = "variant_id", by.y = "variant_id1",
  all.x = FALSE
)
joined_flip <- merge(
  all_snp_dt[!variant_id %in% joined$variant_id],
  eqtl_map[, .(GeneSymbol, ENSG, variant_id1, variant_id2)],
  by.x = "variant_id", by.y = "variant_id2",
  all.x = FALSE
)
all_joined <- rbindlist(list(joined, joined_flip), fill = TRUE)
cat("Variants with gene annotation:", nrow(all_joined), "\n")

# ---------------------------------------------------------------------------
# 6. Load single-ancestry SuSiE credible sets for comparison (EUR-only)
# ---------------------------------------------------------------------------
eur_only_cs <- NULL
if (file.exists(COMBINED_FM)) {
  fm <- fread(COMBINED_FM)
  # Keep EUR SuSiE CS assignments
  if (all(c("chromosome", "position", "allele1", "allele2",
            "susie_cs", "study") %in% names(fm))) {
    eur_only <- fm[grepl("EUR$", study) & susie_cs == 1]
    eur_only[, variant_id := paste(chromosome, position, allele1, allele2, sep = ":")]
    # CS size per locus-study: number of variants per CS group
    eur_only_cs <- eur_only[, .(
      eur_cs_size = .N
    ), by = .(chromosome, locus, study)]
    cat("EUR-only SuSiE loci with CS:", nrow(eur_only_cs), "\n")
  }
}

# ---------------------------------------------------------------------------
# 7. Compute per-gene SuSiEX summary
# ---------------------------------------------------------------------------

# Credible set size per locus (from .cs files)
if (!is.null(all_cs)) {
  all_cs_dt <- as.data.table(all_cs)
  cs_sizes <- all_cs_dt[, .(
    susiex_cs_size_joint = uniqueN(SNP),
    susiex_n_cs          = uniqueN(CS_ID)
  ), by = locus_id]
} else {
  cs_sizes <- data.table(locus_id = character(), susiex_cs_size_joint = integer(),
                          susiex_n_cs = integer())
}

# Gene-level aggregation
gene_summary <- all_joined[, .(
  susiex_max_pip    = max(OVRL_PIP, na.rm = TRUE),
  susiex_n_loci     = uniqueN(locus_id),
  susiex_n_variants = .N,
  trait_pairs       = paste(sort(unique(trait_pair)), collapse = ";"),
  locus_ids         = paste(sort(unique(locus_id)), collapse = ";")
), by = .(GeneSymbol, ENSG)]

# Add CS size per gene (via locus_id join)
gene_cs <- all_joined[OVRL_PIP > 0.1]  # variants with meaningful PIP
gene_cs <- merge(gene_cs, cs_sizes, by = "locus_id", all.x = TRUE)
gene_cs_agg <- gene_cs[, .(
  susiex_n_cs           = max(susiex_n_cs, na.rm = TRUE),
  susiex_cs_size_joint  = min(susiex_cs_size_joint, na.rm = TRUE)  # smallest (tightest) CS
), by = .(GeneSymbol, ENSG)]
gene_cs_agg[is.infinite(susiex_n_cs),          susiex_n_cs := NA_integer_]
gene_cs_agg[is.infinite(susiex_cs_size_joint), susiex_cs_size_joint := NA_integer_]

gene_summary <- merge(gene_summary, gene_cs_agg, by = c("GeneSymbol", "ENSG"), all.x = TRUE)

# Add EUR-only CS size comparison (matched by locus chromosome/position)
# This provides the "power gain" metric: joint CS is tighter than EUR-only
if (!is.null(eur_only_cs)) {
  # Match via shared_loci: find EUR study for each joint locus
  locus_meta <- shared_loci[, c("locus_id", "chr", "window_start", "window_end",
                                 "eur_gwas", "trait_pair")]
  gene_locus <- all_joined[, .(locus_id, GeneSymbol, ENSG, variant_chr = chr, variant_pos = pos)]
  gene_locus <- merge(gene_locus, locus_meta, by = "locus_id")

  # Find EUR-only CS size for same chromosome/locus region
  eur_match <- eur_only_cs[, .(chromosome, locus, study, eur_cs_size)]
  # Match: EUR locus column format is "chr.pos" (e.g. "1.150309638")
  # Extract chr from locus
  eur_match[, eur_chr := as.integer(sub("\\..*", "", locus))]
  eur_match[, eur_pos := as.numeric(sub(".*\\.", "", locus))]

  gene_eur <- merge(
    gene_locus,
    eur_match,
    by.x = c("variant_chr"), by.y = c("eur_chr"),
    allow.cartesian = TRUE
  )
  gene_eur <- gene_eur[abs(eur_pos - variant_pos) < 1e6]  # same region
  eur_cs_per_gene <- gene_eur[, .(
    susiex_cs_size_eur_only = min(eur_cs_size, na.rm = TRUE)
  ), by = .(GeneSymbol, ENSG)]
  eur_cs_per_gene[is.infinite(susiex_cs_size_eur_only), susiex_cs_size_eur_only := NA_integer_]

  gene_summary <- merge(gene_summary, eur_cs_per_gene,
                        by = c("GeneSymbol", "ENSG"), all.x = TRUE)
} else {
  gene_summary[, susiex_cs_size_eur_only := NA_integer_]
}

# Compute CS size reduction
gene_summary[, susiex_cs_reduction := susiex_cs_size_eur_only - susiex_cs_size_joint]

setorder(gene_summary, -susiex_max_pip)

cat("\n=== Top genes by max SuSiEX PIP ===\n")
print(head(gene_summary[, .(GeneSymbol, susiex_max_pip, susiex_n_cs,
                             susiex_cs_size_joint, susiex_cs_size_eur_only,
                             susiex_cs_reduction, trait_pairs)], 20))

# ---------------------------------------------------------------------------
# 8. Write outputs
# ---------------------------------------------------------------------------
out_gene    <- file.path(SUSIEX_DIR, "susiex_gene_summary.csv")
out_variant <- file.path(SUSIEX_DIR, "susiex_variant_summary.csv")

write_csv(as.data.frame(gene_summary), out_gene)
write_csv(as.data.frame(all_snp_dt),   out_variant)

cat("\nOutputs written:\n")
cat("  Gene-level:    ", out_gene,    "\n")
cat("  Variant-level: ", out_variant, "\n")
cat("  Genes with SuSiEX evidence:", nrow(gene_summary), "\n")
cat("  Genes in CS (PIP > 0.5):   ",
    nrow(gene_summary[susiex_max_pip > 0.5]), "\n")
