#!/usr/bin/env Rscript
# 35u_hyprcoloc_broadaway.R
# ---------------------------------------------------------------------------
# Multi-trait colocalization: Broadaway Liver eQTL x ALL GWAS simultaneously
#
# Unlike pairwise COLOC (Scripts 35b/35g/35h), HyPrColoc tests all traits
# simultaneously at each eGene locus. This identifies which GWAS traits share
# a causal variant with the eQTL — far more powerful than pairwise testing.
#
# Algorithm per gene:
#   1. Define cis-window from Broadaway marginal data (variant range + 100kb)
#   2. Extract all variants in cis-window from Broadaway + all GWAS
#   3. Find common SNPs across ALL datasets (by chr:pos matching in hg19)
#   4. Build beta matrix (rows=SNPs, cols=traits) and SE matrix
#   5. Run hyprcoloc::hyprcoloc(betas, ses, trait.names, snp.id)
#   6. Report posterior, trait clusters, regional association probability
#
# Inputs:
#   - data/broadaway_eqtl/Liver_eQTL_Meta_Leads_ST3_20240530.tsv
#   - data/broadaway_eqtl/chr{1-22}_marginal_summary_results.tsv (hg19)
#   - GWAS/MR_Data/hg19/*.tsv.gz (Phase 0D liftover products)
#   - GWAS/MR_Data/GCST90267352_PDFF_Pazoki2022.tsv.gz (natively hg19)
#
# Outputs:
#   - RNA-seq/results/causal_inference/hyprcoloc/hyprcoloc_results.csv
#   - RNA-seq/results/causal_inference/hyprcoloc/hyprcoloc_summary.csv
# ---------------------------------------------------------------------------

.libPaths(c("RNA-seq/.Rlib_causal", .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
})

# Try loading hyprcoloc; if unavailable, use pairwise coloc.abf fallback
HAS_HYPRCOLOC <- requireNamespace("hyprcoloc", quietly = TRUE)
if (HAS_HYPRCOLOC) {
  library(hyprcoloc)
  cat("Using hyprcoloc package for multi-trait colocalization\n")
} else {
  cat("hyprcoloc not available — using pairwise coloc.abf + cluster aggregation\n")
}

cat("=== Script 35u: HyPrColoc Multi-Trait Colocalization ===\n")
cat("Start time:", format(Sys.time()), "\n\n")

# ==============================================================================
# Configuration
# ==============================================================================
BASE_DIR    <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
EQTL_DIR   <- file.path(BASE_DIR, "data/broadaway_eqtl")
GWAS_DIR    <- file.path(BASE_DIR, "GWAS/MR_Data")
HG19_DIR    <- file.path(GWAS_DIR, "hg19")
RESULTS_DIR <- file.path(BASE_DIR, "RNA-seq/results/causal_inference/hyprcoloc")

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

# Broadaway eQTL parameters
EQTL_N <- 1183L

# Minimum shared SNPs per gene to attempt HyPrColoc (raised from 50; see threshold audit)
MIN_SHARED_SNPS <- 100L

# Cis-window padding (bp) beyond observed eQTL variant range
CIS_PADDING <- 100000L

# Checkpoint interval
CHECKPOINT_INTERVAL <- 1000L

# GWAS registry: each entry specifies the file, expected columns, and metadata.
# Coordinate system: "hg19" means positions are in hg19 (matches Broadaway).
# Files in hg19/ are Phase 0D liftover products with hg19_chr and hg19_pos columns.
# PDFF (Pazoki 2022) is natively hg19.
gwas_registry <- list(
  Ghodsian_NAFLD = list(
    file     = file.path(HG19_DIR, "Ghodsian_2021_NAFLD_harmonised_hg19.tsv.gz"),
    format   = "hg19_liftover",
    N        = 778614L,
    type     = "cc",
    s        = 8434 / 778614,  # case proportion
    label    = "Ghodsian NAFLD"
  ),
  UKBB_ALT = list(
    file     = file.path(HG19_DIR, "GCST90019492_UKBB_ALT_harmonised_hg19.tsv.gz"),
    format   = "hg19_liftover",
    N        = 343850L,
    type     = "quant",
    s        = NULL,
    label    = "UKBB ALT"
  ),
  UKBB_AST = list(
    file     = file.path(HG19_DIR, "GCST90019497_UKBB_AST_harmonised_hg19.tsv.gz"),
    format   = "hg19_liftover",
    N        = 343850L,
    type     = "quant",
    s        = NULL,
    label    = "UKBB AST"
  ),
  UKBB_GGT = list(
    file     = file.path(HG19_DIR, "GCST90019507_UKBB_GGT_harmonised_hg19.tsv.gz"),
    format   = "hg19_liftover",
    N        = 343850L,
    type     = "quant",
    s        = NULL,
    label    = "UKBB GGT"
  ),
  PDFF = list(
    file     = file.path(GWAS_DIR, "GCST90267352_PDFF_Pazoki2022.tsv.gz"),
    format   = "pdff_native_hg19",
    N        = 33588L,
    type     = "quant",
    s        = NULL,
    label    = "PDFF (Pazoki 2022)"
  )
)

# ==============================================================================
# 1. Load and validate GWAS datasets
# ==============================================================================
cat("--- Step 1: Loading GWAS summary statistics ---\n")

# Standardized GWAS loader: returns data.table with columns
#   chr (integer), pos_hg19 (integer), effect_allele, other_allele, beta, se
load_gwas <- function(cfg, name) {
  if (!file.exists(cfg$file)) {
    cat("  WARNING:", name, "- file not found:", cfg$file, "\n")
    cat("    Skipping this GWAS.\n")
    return(NULL)
  }

  cat("  Loading", name, "...\n")
  dt <- fread(cfg$file)
  cat("    Raw rows:", nrow(dt), "\n")

  if (cfg$format == "hg19_liftover") {
    # Phase 0D liftover products: original columns + hg19_chr, hg19_pos
    if (!all(c("hg19_chr", "hg19_pos") %in% names(dt))) {
      cat("    WARNING: missing hg19_chr/hg19_pos columns in", name, "\n")
      cat("    Run Phase 0D liftover first. Skipping.\n")
      return(NULL)
    }
    dt[, chr := as.integer(hg19_chr)]
    dt[, pos_hg19 := as.integer(hg19_pos)]

    # Detect beta/SE columns
    beta_col <- intersect(c("beta", "BETA", "hm_beta"), names(dt))[1]
    se_col   <- intersect(c("standard_error", "se", "SE", "sebeta"), names(dt))[1]
    ea_col   <- intersect(c("effect_allele", "alt", "EA", "A1"), names(dt))[1]
    oa_col   <- intersect(c("other_allele", "ref", "NEA", "A2"), names(dt))[1]

    if (is.na(beta_col) || is.na(se_col) || is.na(ea_col) || is.na(oa_col)) {
      cat("    WARNING: cannot detect required columns for", name, ". Skipping.\n")
      return(NULL)
    }

    setnames(dt, c(beta_col, se_col, ea_col, oa_col),
             c("beta", "se", "effect_allele", "other_allele"),
             skip_absent = TRUE)

  } else if (cfg$format == "pdff_native_hg19") {
    # PDFF is natively hg19: chromosome, base_pair_location, ...
    chr_col <- intersect(c("chromosome", "chr", "CHR"), names(dt))[1]
    pos_col <- intersect(c("base_pair_location", "pos", "POS", "BP"), names(dt))[1]
    beta_col <- intersect(c("beta", "BETA"), names(dt))[1]
    se_col   <- intersect(c("standard_error", "se", "SE"), names(dt))[1]
    ea_col   <- intersect(c("effect_allele", "alt", "EA"), names(dt))[1]
    oa_col   <- intersect(c("other_allele", "ref", "NEA"), names(dt))[1]

    if (is.na(chr_col) || is.na(pos_col) || is.na(beta_col) || is.na(se_col) ||
        is.na(ea_col) || is.na(oa_col)) {
      cat("    WARNING: cannot detect required columns for", name, ". Skipping.\n")
      return(NULL)
    }

    dt[, chr := as.integer(gsub("chr", "", get(chr_col), ignore.case = TRUE))]
    dt[, pos_hg19 := as.integer(get(pos_col))]
    setnames(dt, c(beta_col, se_col, ea_col, oa_col),
             c("beta", "se", "effect_allele", "other_allele"),
             skip_absent = TRUE)

  } else {
    cat("    WARNING: unknown format for", name, ". Skipping.\n")
    return(NULL)
  }

  # Filter: valid autosomal, non-missing beta/SE
  dt <- dt[!is.na(chr) & chr %in% 1:22 & !is.na(pos_hg19)]
  dt <- dt[!is.na(beta) & !is.na(se) & se > 0]

  # Detect EAF column for coloc MAF requirement
  eaf_col <- intersect(c("effect_allele_frequency", "EAF", "eaf", "maf", "MAF"), names(dt))[1]
  if (!is.na(eaf_col)) {
    setnames(dt, eaf_col, "effect_allele_frequency", skip_absent = TRUE)
  } else {
    dt[, effect_allele_frequency := NA_real_]
  }

  # Keep only needed columns, deduplicate by chr:pos (lowest |p| proxy: smallest SE)
  dt <- dt[, .(chr, pos_hg19, effect_allele, other_allele, beta, se, effect_allele_frequency)]
  dt[, merge_key := paste0(chr, ":", pos_hg19)]
  dt <- dt[order(merge_key, se)]
  dt <- dt[!duplicated(merge_key)]

  cat("    After QC + dedup:", nrow(dt), "variants\n")
  return(dt)
}

# Load all GWAS
gwas_data <- list()
gwas_meta <- list()

for (gname in names(gwas_registry)) {
  cfg <- gwas_registry[[gname]]
  dt <- load_gwas(cfg, gname)
  if (!is.null(dt)) {
    gwas_data[[gname]] <- dt
    gwas_meta[[gname]] <- cfg
  }
}

n_gwas <- length(gwas_data)
cat("\n  Successfully loaded GWAS:", n_gwas, "/", length(gwas_registry), "\n")

if (n_gwas == 0) {
  cat("ERROR: No GWAS datasets loaded. Cannot proceed.\n")
  cat("  Ensure Phase 0D hg19 liftover files exist in:", HG19_DIR, "\n")
  cat("  Or that PDFF file exists at:", gwas_registry$PDFF$file, "\n")
  quit(save = "no", status = 1)
}

gwas_names <- names(gwas_data)
cat("  Active GWAS:", paste(gwas_names, collapse = ", "), "\n")

# Total traits = eQTL + loaded GWAS
trait_names <- c("Broadaway_eQTL", gwas_names)
n_traits <- length(trait_names)
cat("  Total traits for HyPrColoc:", n_traits, "\n\n")

if (n_traits < 3) {
  cat("WARNING: Only", n_traits, "traits available. HyPrColoc requires >= 2 traits\n")
  cat("  (eQTL + at least 1 GWAS). Proceeding, but multi-trait advantage is limited.\n\n")
}

# ==============================================================================
# 2. Load Broadaway eGene leads + build gene metadata
# ==============================================================================
cat("--- Step 2: Loading Broadaway eGene leads ---\n")

leads <- fread(file.path(EQTL_DIR, "Liver_eQTL_Meta_Leads_ST3_20240530.tsv"))
cat("  Total eQTL signals:", nrow(leads), "\n")

all_eGenes <- unique(leads$Gene)
cat("  Unique eGenes:", length(all_eGenes), "\n")

# Build gene -> chr mapping from leads Variant field (format: CHR_POS_REF_ALT)
leads[, eqtl_chr := as.integer(sub("_.*", "", Variant))]
gene_chr_map <- unique(leads[!is.na(eqtl_chr), .(Gene, chr = eqtl_chr)])
gene_chr_map <- gene_chr_map[, .N, by = .(Gene, chr)][order(-N)][!duplicated(Gene)]
cat("  Gene-to-chromosome mapping:", nrow(gene_chr_map), "genes\n")

# Build gene symbol <-> Ensembl lookup
leads_lookup <- unique(leads[, .(Gene, Ensembl)])

# ==============================================================================
# 3. Load Broadaway marginal eQTL data per chromosome (cached)
# ==============================================================================
cat("\n--- Step 3: Loading Broadaway marginal eQTL data ---\n")

chr_eqtl_cache <- list()

load_chr_eqtl <- function(chr_num) {
  key <- as.character(chr_num)
  if (!is.null(chr_eqtl_cache[[key]])) return(chr_eqtl_cache[[key]])

  fname <- file.path(EQTL_DIR, paste0("chr", chr_num, "_marginal_summary_results.tsv"))
  if (!file.exists(fname)) return(NULL)

  cat("    Loading chr", chr_num, "eQTL marginal data...\n")
  dt <- fread(fname)

  # Broadaway marginal files are hg19 natively: CHR, POS columns
  dt <- dt[!is.na(CHR) & !is.na(POS) & !is.na(Beta) & !is.na(SE) & SE > 0]
  dt[, merge_key := paste0(CHR, ":", POS)]

  # Deduplicate per gene (keep one entry per gene:position)
  dt <- dt[order(GeneSymbol, merge_key, PVAL)]
  dt <- dt[!duplicated(paste(GeneSymbol, merge_key))]

  chr_eqtl_cache[[key]] <<- dt
  cat("      Loaded:", nrow(dt), "variant-gene pairs\n")
  return(dt)
}

# Pre-load all chromosomes
for (chr_num in 1:22) {
  load_chr_eqtl(chr_num)
}

# ==============================================================================
# 4. Allele harmonization
# ==============================================================================
# Harmonize a GWAS dataset's alleles to match the eQTL reference alleles.
# Returns merged data.table with aligned beta values, or NULL if insufficient overlap.

harmonize_to_eqtl <- function(eqtl_dt, gwas_dt) {
  # Merge on chr:pos (hg19)
  # Include EAF from eQTL for coloc MAF requirement
  eqtl_cols <- c("merge_key", "eqtl_EA", "eqtl_NEA", "eqtl_beta", "eqtl_se", "eqtl_eaf")
  eqtl_sub <- eqtl_dt[, .(merge_key, eqtl_EA = EA, eqtl_NEA = NEA,
                            eqtl_beta = Beta, eqtl_se = SE,
                            eqtl_eaf = as.numeric(EAF))]
  gwas_sub <- gwas_dt[, .(merge_key, gwas_EA = effect_allele, gwas_OA = other_allele,
                           gwas_beta = beta, gwas_se = se,
                           gwas_eaf = as.numeric(effect_allele_frequency))]
  merged <- merge(eqtl_sub, gwas_sub, by = "merge_key")

  if (nrow(merged) == 0) return(NULL)

  comp <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C")

  merged[, match_type := "none"]

  # Direct: eQTL EA == GWAS EA and eQTL NEA == GWAS OA
  merged[eqtl_EA == gwas_EA & eqtl_NEA == gwas_OA, match_type := "direct"]

  # Flipped: eQTL EA == GWAS OA and eQTL NEA == GWAS EA
  merged[eqtl_EA == gwas_OA & eqtl_NEA == gwas_EA, match_type := "flipped"]

  # Remove strand-ambiguous unmatched
  merged[, is_ambiguous := (eqtl_EA %in% c("A", "T") & eqtl_NEA %in% c("A", "T")) |
                           (eqtl_EA %in% c("C", "G") & eqtl_NEA %in% c("C", "G"))]
  merged <- merged[!(is_ambiguous & match_type == "none")]

  # Complement matching for remaining
  merged[match_type == "none" & nchar(eqtl_EA) == 1 & nchar(eqtl_NEA) == 1,
    match_type := fifelse(
      comp[eqtl_EA] == gwas_EA & comp[eqtl_NEA] == gwas_OA, "direct",
      fifelse(comp[eqtl_EA] == gwas_OA & comp[eqtl_NEA] == gwas_EA, "flipped",
              "none")
    )]

  # Drop unmatched
  merged <- merged[match_type != "none"]

  # Flip GWAS beta and EAF where alleles are swapped relative to eQTL
  merged[match_type == "flipped", gwas_beta := -gwas_beta]
  merged[match_type == "flipped" & !is.na(gwas_eaf), gwas_eaf := 1 - gwas_eaf]

  merged[, c("is_ambiguous", "match_type") := NULL]
  return(merged)
}

# ==============================================================================
# 5. Run HyPrColoc per eGene
# ==============================================================================
cat("\n--- Step 5: Running HyPrColoc per eGene ---\n")
cat("  eGenes to test:", length(all_eGenes), "\n")
cat("  Traits:", paste(trait_names, collapse = ", "), "\n")
cat("  Min shared SNPs:", MIN_SHARED_SNPS, "\n\n")

results_list <- list()
n_tested     <- 0L
n_skipped_chr <- 0L
n_skipped_eqtl <- 0L
n_skipped_snps <- 0L
n_errors     <- 0L

for (gi in seq_along(all_eGenes)) {
  gene_name <- all_eGenes[gi]

  # Get chromosome
  chr_info <- gene_chr_map[Gene == gene_name]
  if (nrow(chr_info) == 0) {
    n_skipped_chr <- n_skipped_chr + 1L
    next
  }
  gene_chr <- chr_info$chr[1]

  # Get eQTL data for this gene
  eqtl_chr <- chr_eqtl_cache[[as.character(gene_chr)]]
  if (is.null(eqtl_chr)) {
    n_skipped_chr <- n_skipped_chr + 1L
    next
  }

  gene_eqtl <- eqtl_chr[GeneSymbol == gene_name]
  if (nrow(gene_eqtl) == 0) {
    # Try matching on Ensembl ID
    gene_ensg <- leads_lookup$Ensembl[leads_lookup$Gene == gene_name]
    if (length(gene_ensg) > 0) {
      gene_eqtl <- eqtl_chr[ENSG %in% gene_ensg]
    }
  }

  if (nrow(gene_eqtl) < MIN_SHARED_SNPS) {
    n_skipped_eqtl <- n_skipped_eqtl + 1L
    next
  }

  # Define cis-window from eQTL variant range
  pos_range <- range(gene_eqtl$POS)
  cis_start <- max(1L, pos_range[1] - CIS_PADDING)
  cis_end   <- pos_range[2] + CIS_PADDING

  # For each GWAS, harmonize alleles against eQTL within cis-window
  # Collect: list of per-GWAS harmonized tables aligned to eQTL merge_keys
  gwas_harmonized <- list()
  gwas_trait_names <- c()

  for (gname in gwas_names) {
    gdt <- gwas_data[[gname]]

    # Subset to cis-window on correct chromosome
    gdt_cis <- gdt[chr == gene_chr & pos_hg19 >= cis_start & pos_hg19 <= cis_end]

    if (nrow(gdt_cis) < MIN_SHARED_SNPS) next

    harm <- harmonize_to_eqtl(gene_eqtl, gdt_cis)
    if (is.null(harm) || nrow(harm) < MIN_SHARED_SNPS) next

    gwas_harmonized[[gname]] <- harm
    gwas_trait_names <- c(gwas_trait_names, gname)
  }

  if (length(gwas_harmonized) == 0) {
    n_skipped_snps <- n_skipped_snps + 1L
    next
  }

  # Find SNPs common across eQTL AND all harmonized GWAS
  common_snps <- gwas_harmonized[[1]]$merge_key
  if (length(gwas_harmonized) > 1) {
    for (k in 2:length(gwas_harmonized)) {
      common_snps <- intersect(common_snps, gwas_harmonized[[k]]$merge_key)
    }
  }

  if (length(common_snps) < MIN_SHARED_SNPS) {
    n_skipped_snps <- n_skipped_snps + 1L
    next
  }

  # Build beta and SE matrices: rows = SNPs, cols = traits
  # Column order: Broadaway_eQTL, then each GWAS
  active_traits <- c("Broadaway_eQTL", gwas_trait_names)
  n_active <- length(active_traits)

  # eQTL values for common SNPs
  eqtl_sub <- gene_eqtl[merge_key %in% common_snps]
  eqtl_sub <- eqtl_sub[!duplicated(merge_key)]
  setkey(eqtl_sub, merge_key)

  # Enforce consistent SNP ordering
  snp_order <- sort(common_snps)
  eqtl_sub <- eqtl_sub[snp_order]

  beta_mat <- matrix(NA_real_, nrow = length(snp_order), ncol = n_active)
  se_mat   <- matrix(NA_real_, nrow = length(snp_order), ncol = n_active)
  eaf_mat  <- matrix(NA_real_, nrow = length(snp_order), ncol = n_active)
  colnames(beta_mat) <- active_traits
  colnames(se_mat)   <- active_traits
  colnames(eaf_mat)  <- active_traits
  rownames(beta_mat) <- snp_order
  rownames(se_mat)   <- snp_order
  rownames(eaf_mat)  <- snp_order

  # Fill eQTL column
  beta_mat[, 1] <- eqtl_sub$Beta
  se_mat[, 1]   <- eqtl_sub$SE
  eaf_mat[, 1]  <- as.numeric(eqtl_sub$EAF)

  # Fill GWAS columns
  for (k in seq_along(gwas_trait_names)) {
    gname <- gwas_trait_names[k]
    harm <- gwas_harmonized[[gname]]
    harm <- harm[merge_key %in% snp_order]
    harm <- harm[!duplicated(merge_key)]
    setkey(harm, merge_key)
    harm <- harm[snp_order]

    beta_mat[, k + 1] <- harm$gwas_beta
    se_mat[, k + 1]   <- harm$gwas_se
    eaf_mat[, k + 1]  <- harm$gwas_eaf
  }

  # Drop any rows with NA (should be rare after intersection)
  valid_rows <- complete.cases(beta_mat) & complete.cases(se_mat)
  beta_mat <- beta_mat[valid_rows, , drop = FALSE]
  se_mat   <- se_mat[valid_rows, , drop = FALSE]
  eaf_mat  <- eaf_mat[valid_rows, , drop = FALSE]  # FIX (review A08#1): keep eaf_mat aligned to beta/se after NA filter

  if (nrow(beta_mat) < MIN_SHARED_SNPS) {
    n_skipped_snps <- n_skipped_snps + 1L
    next
  }

  # Run multi-trait colocalization
  if (HAS_HYPRCOLOC) {
    # Use hyprcoloc package directly
    res <- tryCatch({
      hyprcoloc(
        effect.est  = beta_mat,
        effect.se   = se_mat,
        trait.names = active_traits,
        snp.id      = rownames(beta_mat)
      )
    }, error = function(e) {
      if (n_errors < 10) cat("  ERROR (", gene_name, "):", conditionMessage(e), "\n")
      n_errors <<- n_errors + 1L
      NULL
    })

    if (is.null(res)) next
    n_tested <- n_tested + 1L
    coloc_method <- "hyprcoloc"  # FIX (review A08#2): real HyPrColoc cluster posterior

    hypr_res <- res$results
    if (nrow(hypr_res) == 0 || all(is.na(hypr_res$traits))) {
      trait_str <- NA_character_; n_cluster_traits <- 0L
      pp <- 0; rp <- 0; csnp <- NA_character_; csnp_prob <- 0
    } else {
      top <- hypr_res[1, , drop = FALSE]
      trait_str <- as.character(top$traits)
      if (is.na(trait_str) || trait_str == "" || trait_str == "None") {
        n_cluster_traits <- 0L; trait_str <- NA_character_
      } else {
        n_cluster_traits <- length(strsplit(trait_str, ", ")[[1]])
      }
      pp <- as.numeric(top$posterior_prob)
      rp <- as.numeric(top$regional_prob)
      csnp <- as.character(top$candidate_snp)
      csnp_prob <- as.numeric(top$posterior_explained_by_snp)
    }

  } else {
    # Fallback: pairwise coloc.abf eQTL vs each GWAS, then aggregate
    # Multi-trait posterior = product of pairwise PP.H4 (conservative)
    # Trait cluster = all GWAS with PP.H4 > 0.5 against eQTL

    pairwise_pp4 <- numeric(0)
    coloc_traits  <- character(0)
    best_snp      <- NA_character_
    best_pp4      <- 0

    for (k in seq_along(gwas_trait_names)) {
      gname <- gwas_trait_names[k]
      gcfg  <- gwas_meta[[gname]]

      # coloc.abf needs: beta, varbeta, snp, position, type, N, MAF (or sdY)
      eqtl_maf <- eaf_mat[, 1]
      eqtl_maf <- pmin(eqtl_maf, 1 - eqtl_maf)  # convert EAF to MAF
      eqtl_list <- list(
        beta    = beta_mat[, 1],
        varbeta = se_mat[, 1]^2,
        snp     = rownames(beta_mat),
        position = as.integer(sub(".*:", "", rownames(beta_mat))),
        type    = "quant",
        N       = EQTL_N,
        MAF     = eqtl_maf
      )
      gwas_list <- list(
        beta    = beta_mat[, k + 1],
        varbeta = se_mat[, k + 1]^2,
        snp     = rownames(beta_mat),
        position = as.integer(sub(".*:", "", rownames(beta_mat))),
        type    = gcfg$type,
        N       = gcfg$N
      )
      if (gcfg$type == "quant") gwas_list$sdY <- 1
      if (gcfg$type == "cc" && !is.null(gcfg$s)) gwas_list$s <- gcfg$s

      coloc_res <- tryCatch(coloc.abf(eqtl_list, gwas_list, p12 = 1e-5),
                            error = function(e) {
                              if (n_errors < 5) cat("  COLOC ERROR (", gene_name, "/", gname, "):", conditionMessage(e), "\n")
                              NULL
                            })
      if (is.null(coloc_res)) next

      pp4 <- coloc_res$summary["PP.H4.abf"]
      pairwise_pp4[gname] <- pp4

      if (pp4 > 0.5) coloc_traits <- c(coloc_traits, gname)
      if (pp4 > best_pp4) {
        best_pp4 <- pp4
        # Get top SNP from coloc results
        if (!is.null(coloc_res$results)) {
          top_row <- which.max(coloc_res$results$SNP.PP.H4)
          best_snp <- coloc_res$results$snp[top_row]
        }
      }
    }

    if (length(pairwise_pp4) == 0) {
      n_errors <- n_errors + 1L
      next
    }
    n_tested <- n_tested + 1L
    coloc_method <- "coloc_abf_gmean_fallback"  # FIX (review A08#2): NOT a HyPrColoc posterior — ad-hoc geometric mean of pairwise ABF PP.H4

    # Multi-trait posterior: geometric mean of pairwise PP.H4 for colocalizing traits
    if (length(coloc_traits) > 0) {
      pp <- exp(mean(log(pairwise_pp4[coloc_traits])))
      trait_str <- paste(c("Broadaway_eQTL", coloc_traits), collapse = ", ")
      n_cluster_traits <- length(coloc_traits) + 1L  # +1 for eQTL
    } else {
      pp <- max(pairwise_pp4, na.rm = TRUE)
      trait_str <- NA_character_
      n_cluster_traits <- 0L
    }
    rp <- pp  # approximate
    csnp <- best_snp
    csnp_prob <- best_pp4
  }

  results_list[[gene_name]] <- data.table(
    gene                        = gene_name,
    gene_symbol                 = gene_name,
    ensembl                     = paste(unique(gene_eqtl$ENSG[gene_eqtl$ENSG != ""]),
                                        collapse = ";"),
    chr                         = gene_chr,
    n_snps                      = nrow(beta_mat),
    n_traits_tested             = n_active,
    traits_tested               = paste(active_traits, collapse = ";"),
    hyprcoloc_posterior         = pp,
    hyprcoloc_regional_prob     = rp,
    hyprcoloc_n_traits          = n_cluster_traits,
    hyprcoloc_traits_list       = trait_str,
    hyprcoloc_candidate_snp     = csnp,
    hyprcoloc_candidate_snp_prob = csnp_prob,
    method                      = coloc_method  # FIX (review A08#2): hyprcoloc vs coloc_abf_gmean_fallback
  )

  # Progress
  if (gi %% 500 == 0) {
    cat("  [", gi, "/", length(all_eGenes), "] tested:", n_tested,
        " skipped_chr:", n_skipped_chr,
        " skipped_eqtl:", n_skipped_eqtl,
        " skipped_snps:", n_skipped_snps,
        " errors:", n_errors, "\n")
  }

  # Checkpoint: save intermediate results
  if (n_tested > 0 && n_tested %% CHECKPOINT_INTERVAL == 0) {
    cat("  >>> Checkpoint at", n_tested, "tested genes\n")
    checkpoint_dt <- rbindlist(results_list, fill = TRUE)
    fwrite(checkpoint_dt, file.path(RESULTS_DIR, "hyprcoloc_results_checkpoint.csv"))
  }
}

# ==============================================================================
# 6. Compile and save final results
# ==============================================================================
cat("\n--- Step 6: Compiling results ---\n")

cat("  Genes tested:", n_tested, "\n")
cat("  Genes skipped (no chr / no eQTL data):", n_skipped_chr, "\n")
cat("  Genes skipped (< ", MIN_SHARED_SNPS, " eQTL variants):", n_skipped_eqtl, "\n")
cat("  Genes skipped (< ", MIN_SHARED_SNPS, " shared SNPs across traits):", n_skipped_snps, "\n")
cat("  Genes with errors:", n_errors, "\n")

if (length(results_list) == 0) {
  cat("  No HyPrColoc results. Writing empty output.\n")
  fwrite(data.table(
    gene = character(), gene_symbol = character(), chr = integer(),
    n_snps = integer(), hyprcoloc_posterior = numeric(),
    hyprcoloc_regional_prob = numeric(), hyprcoloc_n_traits = integer(),
    hyprcoloc_traits_list = character(), hyprcoloc_candidate_snp = character(),
    hyprcoloc_candidate_snp_prob = numeric(), method = character()
  ), file.path(RESULTS_DIR, "hyprcoloc_results.csv"))
  quit(save = "no", status = 0)
}

hyprcoloc_res <- rbindlist(results_list, fill = TRUE)
setorder(hyprcoloc_res, -hyprcoloc_posterior)

fwrite(hyprcoloc_res, file.path(RESULTS_DIR, "hyprcoloc_results.csv"))

# Clean up checkpoint file
checkpoint_file <- file.path(RESULTS_DIR, "hyprcoloc_results_checkpoint.csv")
if (file.exists(checkpoint_file)) {
  file.remove(checkpoint_file)
}

# ==============================================================================
# 7. Summary statistics
# ==============================================================================
cat("\n--- Step 7: Summary ---\n\n")

# Overall
cat("  HyPrColoc Summary:\n")
cat("    Traits:", paste(trait_names, collapse = ", "), "\n")
cat("    Genes tested:", n_tested, "\n")

# Colocalization thresholds
n_coloc_any   <- sum(hyprcoloc_res$hyprcoloc_posterior > 0, na.rm = TRUE)
n_coloc_05    <- sum(hyprcoloc_res$hyprcoloc_posterior > 0.5, na.rm = TRUE)
n_coloc_08    <- sum(hyprcoloc_res$hyprcoloc_posterior > 0.8, na.rm = TRUE)
n_coloc_025   <- sum(hyprcoloc_res$hyprcoloc_posterior > 0.25, na.rm = TRUE)

cat("    PP > 0 (any cluster):", n_coloc_any, "\n")
cat("    PP > 0.25:", n_coloc_025, "\n")
cat("    PP > 0.5:", n_coloc_05, "\n")
cat("    PP > 0.8:", n_coloc_08, "\n")

# Trait cluster sizes
coloc_hits <- hyprcoloc_res[hyprcoloc_posterior > 0.25 & !is.na(hyprcoloc_traits_list)]
if (nrow(coloc_hits) > 0) {
  cat("\n  Cluster size distribution (PP > 0.25):\n")
  cluster_tab <- table(coloc_hits$hyprcoloc_n_traits)
  for (sz in sort(as.integer(names(cluster_tab)))) {
    cat("    ", sz, "traits:", cluster_tab[as.character(sz)], "genes\n")
  }

  # Which GWAS appear most often in clusters
  all_cluster_traits <- unlist(strsplit(coloc_hits$hyprcoloc_traits_list, ", "))
  trait_freq <- sort(table(all_cluster_traits), decreasing = TRUE)
  cat("\n  Trait frequency in clusters (PP > 0.25):\n")
  for (tn in names(trait_freq)) {
    cat("    ", tn, ":", trait_freq[tn], "\n")
  }
}

# Top hits
top_hits <- hyprcoloc_res[hyprcoloc_posterior > 0.5]
if (nrow(top_hits) > 0) {
  cat("\n  Top HyPrColoc hits (PP > 0.5):\n")
  for (i in seq_len(min(nrow(top_hits), 30))) {
    cat("    ", top_hits$gene[i],
        ": PP =", round(top_hits$hyprcoloc_posterior[i], 3),
        ", traits =", top_hits$hyprcoloc_n_traits[i],
        " (", top_hits$hyprcoloc_traits_list[i], ")",
        ", SNP =", top_hits$hyprcoloc_candidate_snp[i], "\n")
  }
}

# Save summary
summary_dt <- data.table(
  metric = c("n_gwas_loaded", "n_traits_total", "trait_names",
             "genes_tested", "genes_skipped_chr", "genes_skipped_eqtl",
             "genes_skipped_snps", "genes_errors",
             "min_shared_snps", "cis_padding_bp",
             "PP_gt_0", "PP_gt_0.25", "PP_gt_0.5", "PP_gt_0.8"),
  value = c(n_gwas, n_traits, paste(trait_names, collapse = ";"),
            n_tested, n_skipped_chr, n_skipped_eqtl,
            n_skipped_snps, n_errors,
            MIN_SHARED_SNPS, CIS_PADDING,
            n_coloc_any, n_coloc_025, n_coloc_05, n_coloc_08)
)
fwrite(summary_dt, file.path(RESULTS_DIR, "hyprcoloc_summary.csv"))

# ==============================================================================
# 8. Validation: known MASLD genes
# ==============================================================================
cat("\n--- Step 8: Validation against known genes ---\n")

known_genes <- c("HSD17B13", "PNPLA3", "TM6SF2", "MBOAT7", "GCKR",
                 "MARC1", "THRB", "DGAT2", "FXR", "NR1H4", "PPARA",
                 "SLC39A8", "SORT1", "CELSR2", "RORA", "CDK6")

for (g in known_genes) {
  hit <- hyprcoloc_res[gene == g]
  if (nrow(hit) > 0) {
    cat("  ", g, ": PP =", round(hit$hyprcoloc_posterior[1], 3),
        ", traits =", hit$hyprcoloc_n_traits[1],
        " (", hit$hyprcoloc_traits_list[1], ")",
        ", n_snps =", hit$n_snps[1], "\n")
  } else {
    cat("  ", g, ": not tested (no Broadaway eQTL)\n")
  }
}

# ==============================================================================
# 9. Compare with pairwise COLOC
# ==============================================================================
cat("\n--- Step 9: Comparison with pairwise COLOC ---\n")

pairwise_dirs <- list(
  UKBB_ALT = file.path(BASE_DIR, "RNA-seq/results/causal_inference/broadaway_ukbb"),
  UKBB_AST = file.path(BASE_DIR, "RNA-seq/results/causal_inference/broadaway_ukbb_ast"),
  UKBB_GGT = file.path(BASE_DIR, "RNA-seq/results/causal_inference/broadaway_ukbb_ggt"),
  PDFF     = file.path(BASE_DIR, "RNA-seq/results/causal_inference/broadaway_pdff")
)

for (pname in names(pairwise_dirs)) {
  pfile <- file.path(pairwise_dirs[[pname]], "coloc_results.csv")
  if (!file.exists(pfile)) next

  pres <- fread(pfile)
  if (nrow(pres) == 0 || !"PP.H4" %in% names(pres)) next

  # Genes with pairwise PP.H4 > 0.5
  pairwise_hits <- pres[PP.H4 > 0.5]$gene
  # Genes with HyPrColoc PP > 0.5 and this GWAS in the cluster
  hypr_hits <- hyprcoloc_res[hyprcoloc_posterior > 0.5 &
                               grepl(pname, hyprcoloc_traits_list, fixed = TRUE)]$gene

  shared <- intersect(pairwise_hits, hypr_hits)
  pairwise_only <- setdiff(pairwise_hits, hypr_hits)
  hypr_only <- setdiff(hypr_hits, pairwise_hits)

  cat("  ", pname, ":\n")
  cat("    Pairwise PP.H4 > 0.5:", length(pairwise_hits), "\n")
  cat("    HyPrColoc PP > 0.5 (includes", pname, "):", length(hypr_hits), "\n")
  cat("    Both:", length(shared),
      " | Pairwise-only:", length(pairwise_only),
      " | HyPrColoc-only:", length(hypr_only), "\n")
}

cat("\n=== Script 35u complete ===\n")
cat("End time:", format(Sys.time()), "\n")
cat("Results:", file.path(RESULTS_DIR, "hyprcoloc_results.csv"), "\n")
