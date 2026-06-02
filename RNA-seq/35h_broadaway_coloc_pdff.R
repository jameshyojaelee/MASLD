#!/usr/bin/env Rscript
# 35h_broadaway_coloc_pdff.R
# ---------------------------------------------------------------------------
# Colocalization: Broadaway Liver eQTL (N=1,183) x PDFF GWAS (N=33,588)
#
# Key differences from 35b/35g:
#   - PDFF GWAS (Pazoki 2022, GCST90267352) is GRCh37 — same as Broadaway native
#   - NO liftover needed: merge eQTL and GWAS in native hg19 coordinate space
#   - Lower N (33K vs 344K) limits power; expect fewer but steatosis-specific hits
#   - GWAS-SSF format may differ from UKBB harmonised — handles column mapping
#
# Inputs:
#   - data/broadaway_eqtl/Liver_eQTL_Meta_Leads_ST3_20240530.tsv
#   - data/broadaway_eqtl/chr{1-22}_marginal_summary_results.tsv
#   - GWAS/MR_Data/GCST90267352_PDFF_Pazoki2022.tsv.gz
#
# Outputs:
#   - RNA-seq/results/causal_inference/broadaway_pdff/coloc_results.csv
#   - RNA-seq/results/causal_inference/broadaway_pdff/coloc_summary.csv
#   - RNA-seq/results/causal_inference/broadaway_pdff/coloc_sensitivity/
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(coloc)
  library(data.table)
  library(ggplot2)
})

cat("=== Script 35h: Broadaway eQTL COLOC x PDFF GWAS (Pazoki 2022) ===\n")
cat("Start time:", format(Sys.time()), "\n\n")

# ==============================================================================
# Configuration
# ==============================================================================
BASE_DIR    <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
EQTL_DIR   <- file.path(BASE_DIR, "data/broadaway_eqtl")
GWAS_FILE   <- file.path(BASE_DIR, "GWAS/MR_Data/GCST90267352_PDFF_Pazoki2022.tsv.gz")
RESULTS_DIR <- file.path(BASE_DIR, "RNA-seq/results/causal_inference/broadaway_pdff")
SENS_DIR    <- file.path(RESULTS_DIR, "coloc_sensitivity")

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SENS_DIR, recursive = TRUE, showWarnings = FALSE)

# GWAS parameters (PDFF: quantitative, N=33,588)
GWAS_N    <- 33588L
GWAS_TYPE <- "quant"

# Broadaway eQTL parameters
EQTL_N    <- 1183L

# Minimum overlapping SNPs for COLOC (raised from 10; see threshold audit)
MIN_SNPS <- 100L

# COLOC priors
COLOC_P1  <- 1e-4
COLOC_P2  <- 1e-4
COLOC_P12 <- 5e-6

# ==============================================================================
# 1. Load PDFF GWAS (GRCh37 — NO liftover needed)
# ==============================================================================
cat("--- Step 1: Loading PDFF GWAS (Pazoki 2022, GRCh37) ---\n")

if (!file.exists(GWAS_FILE)) {
  stop("GWAS file not found: ", GWAS_FILE,
       "\nRun GWAS/download_ukbb_liver_enzymes.sh first.")
}

gwas <- fread(GWAS_FILE)
cat("  Raw GWAS rows:", nrow(gwas), "\n")
cat("  Columns:", paste(names(gwas), collapse = ", "), "\n")

# Handle GWAS-SSF format column mapping
# Standard GWAS-SSF: chromosome, base_pair_location, effect_allele, other_allele, beta, standard_error, p_value
# Some variants: chr, pos, ref, alt, etc.
col_map <- list(
  chr_col = NULL, pos_col = NULL, beta_col = NULL, se_col = NULL, p_col = NULL,
  ea_col = NULL, oa_col = NULL
)

# Detect chromosome column
for (c in c("chromosome", "chr", "CHR", "hm_chrom", "#CHROM")) {
  if (c %in% names(gwas)) { col_map$chr_col <- c; break }
}
# Detect position column
for (c in c("base_pair_location", "pos", "POS", "BP", "hm_pos")) {
  if (c %in% names(gwas)) { col_map$pos_col <- c; break }
}
# Detect beta column
for (c in c("beta", "BETA", "hm_beta", "effect_size")) {
  if (c %in% names(gwas)) { col_map$beta_col <- c; break }
}
# Detect SE column
for (c in c("standard_error", "se", "SE", "hm_se")) {
  if (c %in% names(gwas)) { col_map$se_col <- c; break }
}
# Detect p-value column
for (c in c("p_value", "pvalue", "pval", "P", "hm_p")) {
  if (c %in% names(gwas)) { col_map$p_col <- c; break }
}
# Detect effect allele column
for (c in c("effect_allele", "EA", "A1", "hm_effect_allele", "alt")) {
  if (c %in% names(gwas)) { col_map$ea_col <- c; break }
}
# Detect other allele column
for (c in c("other_allele", "NEA", "A2", "hm_other_allele", "ref")) {
  if (c %in% names(gwas)) { col_map$oa_col <- c; break }
}

# Verify all columns found
missing <- names(col_map)[sapply(col_map, is.null)]
if (length(missing) > 0) {
  stop("Could not detect required GWAS columns: ", paste(missing, collapse = ", "),
       "\nAvailable columns: ", paste(names(gwas), collapse = ", "))
}

cat("  Column mapping:\n")
for (nm in names(col_map)) {
  cat("    ", nm, "→", col_map[[nm]], "\n")
}

# Standardize column names
setnames(gwas, col_map$chr_col, "chr_raw")
setnames(gwas, col_map$pos_col, "pos_hg19")
setnames(gwas, col_map$beta_col, "beta")
setnames(gwas, col_map$se_col, "standard_error")
setnames(gwas, col_map$p_col, "p_value")
setnames(gwas, col_map$ea_col, "effect_allele")
setnames(gwas, col_map$oa_col, "other_allele")

# Filter to valid autosomal rows
gwas[, chr := as.integer(gsub("chr", "", chr_raw, ignore.case = TRUE))]
gwas[, pos_hg19 := as.integer(pos_hg19)]
gwas <- gwas[!is.na(chr) & chr %in% 1:22 & !is.na(pos_hg19) & !is.na(p_value)]

# Need beta + SE for COLOC
gwas <- gwas[!is.na(beta) & !is.na(standard_error) & standard_error > 0]

# Deduplicate by position (keep lowest p-value)
gwas <- gwas[order(chr, pos_hg19, p_value)]
gwas <- gwas[!duplicated(paste(chr, pos_hg19))]
gwas[, merge_key := paste0(chr, ":", pos_hg19)]

cat("  GWAS after dedup:", nrow(gwas), "variants\n")
cat("  Genome-wide significant (p<5e-8):", sum(gwas$p_value < 5e-8), "\n")

# ==============================================================================
# 2. Load Broadaway leads — ALL eGenes (no DEG filter)
# ==============================================================================
cat("\n--- Step 2: Loading Broadaway eGene leads ---\n")

leads <- fread(file.path(EQTL_DIR, "Liver_eQTL_Meta_Leads_ST3_20240530.tsv"))
cat("  Total eQTL signals:", nrow(leads), "\n")

all_eGenes <- unique(leads$Gene)
cat("  Unique eGenes (ALL, no DEG filter):", length(all_eGenes), "\n")

# Build gene → chr mapping from leads Variant field (format: CHR_POS_REF_ALT)
leads[, eqtl_chr := as.integer(sub("_.*", "", Variant))]
gene_chr_map <- unique(leads[!is.na(eqtl_chr), .(Gene, chr = eqtl_chr)])
gene_chr_map <- gene_chr_map[, .N, by = .(Gene, chr)][order(-N)][!duplicated(Gene)]
cat("  Gene-to-chromosome mapping:", nrow(gene_chr_map), "genes\n")

# Build gene symbol <-> Ensembl lookup from leads
leads_lookup <- unique(leads[, .(Gene, Ensembl)])

# ==============================================================================
# 3. Helper functions
# ==============================================================================

# Per-chromosome eQTL cache — NO liftover, use native hg19 POS
chr_eqtl_cache <- list()

load_chr_eqtl <- function(chr_num) {
  key <- as.character(chr_num)
  if (!is.null(chr_eqtl_cache[[key]])) return(chr_eqtl_cache[[key]])

  fname <- file.path(EQTL_DIR, paste0("chr", chr_num, "_marginal_summary_results.tsv"))
  if (!file.exists(fname)) return(NULL)

  cat("    Loading chr", chr_num, "eQTL data (hg19, no liftover)...\n")
  dt <- fread(fname)

  # Merge key uses native hg19 coordinates (both PDFF and Broadaway are hg19)
  dt[, merge_key := paste0(CHR, ":", POS)]

  chr_eqtl_cache[[key]] <<- dt
  return(dt)
}

# Allele harmonization
harmonize_alleles <- function(dt) {
  dt <- copy(dt)
  comp <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C")

  dt[, match_type := "none"]

  # Broadaway: EA = effect allele, NEA = non-effect allele
  # PDFF: effect_allele, other_allele
  dt[EA == effect_allele & NEA == other_allele, match_type := "direct"]
  dt[EA == other_allele & NEA == effect_allele, match_type := "flipped"]

  # Remove strand-ambiguous variants that failed direct/flipped match
  # (complement matching unreliable for A/T and C/G pairs)
  dt[, is_ambiguous := (EA %in% c("A","T") & NEA %in% c("A","T")) |
                        (EA %in% c("C","G") & NEA %in% c("C","G"))]
  n_ambig_unmatched <- sum(dt$is_ambiguous & dt$match_type == "none", na.rm = TRUE)
  cat("    Strand-ambiguous unmatched:", n_ambig_unmatched, "\n")
  dt <- dt[!(is_ambiguous & match_type == "none")]

  # Complement matching for remaining non-ambiguous strand SNPs
  dt[match_type == "none" & nchar(EA) == 1 & nchar(NEA) == 1,
     match_type := fifelse(
       comp[EA] == effect_allele & comp[NEA] == other_allele, "direct",
       fifelse(comp[EA] == other_allele & comp[NEA] == effect_allele, "flipped",
               "none")
     )]

  n_direct  <- sum(dt$match_type == "direct")
  n_flipped <- sum(dt$match_type == "flipped")
  n_none    <- sum(dt$match_type == "none")
  cat("    Allele harmonization: direct=", n_direct, " flipped=", n_flipped,
      " dropped=", n_none, "\n")

  dt <- dt[match_type != "none"]

  if (any(dt$match_type == "flipped")) {
    dt[match_type == "flipped", Beta := -Beta]
  }

  dt[, is_ambiguous := NULL]
  return(dt)
}

# ==============================================================================
# 4. Run COLOC per eGene
# ==============================================================================
cat("\n--- Step 4: Running COLOC ---\n")

candidate_chrs <- unique(gene_chr_map[Gene %in% all_eGenes]$chr)
cat("  Chromosomes to process:", paste(sort(candidate_chrs), collapse = ", "), "\n")
cat("  Total eGenes to test:", length(all_eGenes), "\n\n")

results_list <- list()
n_tested     <- 0L
n_skipped    <- 0L
n_low_snps   <- 0L
n_no_harm    <- 0L

for (gene_name in all_eGenes) {
  chr_info <- gene_chr_map[Gene == gene_name]
  if (nrow(chr_info) == 0) {
    n_skipped <- n_skipped + 1L
    next
  }
  gene_chr <- chr_info$chr[1]

  eqtl_data <- load_chr_eqtl(gene_chr)
  if (is.null(eqtl_data)) {
    n_skipped <- n_skipped + 1L
    next
  }

  gene_eqtl <- eqtl_data[GeneSymbol == gene_name]
  if (nrow(gene_eqtl) == 0) {
    gene_ensg <- leads_lookup$Ensembl[leads_lookup$Gene == gene_name]
    if (length(gene_ensg) > 0) {
      gene_eqtl <- eqtl_data[ENSG %in% gene_ensg]
    }
  }

  if (nrow(gene_eqtl) < MIN_SNPS) {
    n_low_snps <- n_low_snps + 1L
    next
  }

  # Merge eQTL with GWAS on chr:pos (hg19 — both datasets native)
  merged <- merge(gene_eqtl, gwas, by = "merge_key", suffixes = c(".eqtl", ".gwas"))

  if (nrow(merged) < MIN_SNPS) {
    n_low_snps <- n_low_snps + 1L
    next
  }

  merged <- merged[!is.na(Beta) & !is.na(SE) & SE > 0 &
                   !is.na(beta) & !is.na(standard_error) & standard_error > 0]

  if (nrow(merged) < MIN_SNPS) {
    n_low_snps <- n_low_snps + 1L
    next
  }

  merged <- harmonize_alleles(merged)

  if (nrow(merged) < MIN_SNPS) {
    n_no_harm <- n_no_harm + 1L
    next
  }

  dataset1 <- list(
    snp     = merged$merge_key,
    beta    = merged$beta,
    varbeta = (merged$standard_error)^2,
    type    = GWAS_TYPE,
    sdY     = 1,
    N       = GWAS_N
  )

  dataset2 <- list(
    snp     = merged$merge_key,
    beta    = merged$Beta,
    varbeta = (merged$SE)^2,
    type    = "quant",
    sdY     = 1,
    N       = EQTL_N
  )

  my.res <- tryCatch(
    coloc.abf(dataset1, dataset2, p1 = COLOC_P1, p2 = COLOC_P2, p12 = COLOC_P12),
    error = function(e) {
      cat("  coloc.abf error:", gene_name, "-", conditionMessage(e), "\n")
      NULL
    }
  )

  if (!is.null(my.res)) {
    n_tested <- n_tested + 1L

    top_snp_row <- my.res$results[which.max(my.res$results$SNP.PP.H4), ]

    results_list[[gene_name]] <- data.table(
      gene        = gene_name,
      ensembl     = paste(unique(gene_eqtl$ENSG[gene_eqtl$ENSG != ""]), collapse = ";"),
      chr         = gene_chr,
      PP.H0       = my.res$summary["PP.H0.abf"],
      PP.H1       = my.res$summary["PP.H1.abf"],
      PP.H2       = my.res$summary["PP.H2.abf"],
      PP.H3       = my.res$summary["PP.H3.abf"],
      PP.H4       = my.res$summary["PP.H4.abf"],
      n_snps      = nrow(merged),
      n_eqtl_snps = nrow(gene_eqtl),
      top_snp     = top_snp_row$snp,
      top_snp_PP  = top_snp_row$SNP.PP.H4,
      eqtl_source = "Broadaway2024",
      eqtl_N      = EQTL_N,
      gwas_source = "PDFF_Pazoki2022",
      gwas_N      = GWAS_N,
      method      = "abf_beta_varbeta"
    )

    if (my.res$summary["PP.H4.abf"] > 0.3) {
      tryCatch({
        pdf(file.path(SENS_DIR, paste0("sensitivity_", gene_name, ".pdf")),
            width = 10, height = 8)
        sensitivity(my.res, rule = "H4 > 0.5")
        dev.off()
      }, error = function(e) {
        tryCatch(dev.off(), error = function(e2) NULL)
      })
    }

    if (n_tested %% 200 == 0) {
      cat("  [", n_tested, "tested,", n_low_snps, "low-SNP,",
          n_no_harm, "no-harmonize,", n_skipped, "skipped]\n")
    }
  }
}

# ==============================================================================
# 5. Compile and save results
# ==============================================================================
cat("\n--- Step 5: Compiling results ---\n")

cat("  Genes tested:", n_tested, "\n")
cat("  Genes skipped (no chr mapping):", n_skipped, "\n")
cat("  Genes skipped (<", MIN_SNPS, "shared SNPs):", n_low_snps, "\n")
cat("  Genes skipped (allele harmonization):", n_no_harm, "\n")

if (length(results_list) == 0) {
  cat("  No COLOC results.\n")
  fwrite(data.table(gene = character(), PP.H4 = numeric()),
         file.path(RESULTS_DIR, "coloc_results.csv"))
  quit(save = "no", status = 0)
}

coloc_res <- rbindlist(results_list, fill = TRUE)
setorder(coloc_res, -PP.H4)

fwrite(coloc_res, file.path(RESULTS_DIR, "coloc_results.csv"))

cat("\n  COLOC Summary (Broadaway N=1,183 x PDFF N=33,588):\n")
cat("    Genes tested:", n_tested, "\n")
cat("    PP.H4 > 0.8:", nrow(coloc_res[PP.H4 > 0.8]), "\n")
cat("    PP.H4 > 0.5:", nrow(coloc_res[PP.H4 > 0.5]), "\n")
cat("    PP.H4 > 0.3:", nrow(coloc_res[PP.H4 > 0.3]), "\n")

top_hits <- coloc_res[PP.H4 > 0.3]
if (nrow(top_hits) > 0) {
  cat("\n  Top COLOC hits (PP.H4 > 0.3):\n")
  for (i in seq_len(min(nrow(top_hits), 30))) {
    cat("    ", top_hits$gene[i], ": PP.H4 =", round(top_hits$PP.H4[i], 3),
        "(", top_hits$n_snps[i], "SNPs, top:", top_hits$top_snp[i], ")\n")
  }
}

summary_dt <- data.table(
  metric = c("genes_tested", "genes_skipped_no_chr", "genes_skipped_low_snps",
             "genes_skipped_harmonization",
             "PP.H4_gt_0.8", "PP.H4_gt_0.5", "PP.H4_gt_0.3",
             "eqtl_source", "eqtl_N", "gwas_source", "gwas_N",
             "coordinate_space", "liftover"),
  value = c(n_tested, n_skipped, n_low_snps, n_no_harm,
            nrow(coloc_res[PP.H4 > 0.8]),
            nrow(coloc_res[PP.H4 > 0.5]),
            nrow(coloc_res[PP.H4 > 0.3]),
            "Broadaway2024", EQTL_N,
            "PDFF_GCST90267352", GWAS_N,
            "hg19", "none_both_native")
)
fwrite(summary_dt, file.path(RESULTS_DIR, "coloc_summary.csv"))

# ==============================================================================
# 6. Validation + comparison
# ==============================================================================
cat("\n--- Step 6: Validation ---\n")

# Known steatosis-associated genes
steatosis_genes <- c("PNPLA3", "TM6SF2", "HSD17B13", "GCKR", "MBOAT7",
                     "MARC1", "SPTLC3", "FABP1", "CIDEB", "FASN")
cat("\n  Steatosis gene validation:\n")
for (g in steatosis_genes) {
  hit <- coloc_res[gene == g]
  if (nrow(hit) > 0) {
    cat("    ", g, ": PP.H4 =", round(hit$PP.H4[1], 3),
        "(", hit$n_snps[1], "SNPs)\n")
  } else {
    cat("    ", g, ": not tested (no eQTL or insufficient overlap)\n")
  }
}

# Compare with ALT results
alt_file <- file.path(BASE_DIR, "RNA-seq/results/causal_inference/broadaway_ukbb/coloc_results.csv")
if (file.exists(alt_file)) {
  alt_res <- fread(alt_file)
  if (nrow(alt_res) > 0 && "PP.H4" %in% names(alt_res)) {
    cat("\n  Comparison with Broadaway x UKBB ALT:\n")
    cat("    ALT PP.H4 > 0.5:", nrow(alt_res[PP.H4 > 0.5]), "\n")
    cat("    PDFF PP.H4 > 0.5:", nrow(coloc_res[PP.H4 > 0.5]), "\n")
    shared <- intersect(alt_res$gene, coloc_res$gene)
    if (length(shared) > 5) {
      comp <- merge(
        alt_res[, .(gene, PP.H4_alt = PP.H4)],
        coloc_res[, .(gene, PP.H4_pdff = PP.H4)],
        by = "gene"
      )
      cat("    Shared genes tested:", nrow(comp), "\n")
      cat("    Pearson r (PP.H4):", round(cor(comp$PP.H4_alt, comp$PP.H4_pdff,
                                               use = "complete.obs"), 3), "\n")
      # PDFF-unique hits (steatosis-specific)
      pdff_unique <- setdiff(coloc_res[PP.H4 > 0.5]$gene, alt_res[PP.H4 > 0.5]$gene)
      cat("    PDFF-unique (PP.H4>0.5):", length(pdff_unique),
          "— steatosis-specific candidates\n")
      if (length(pdff_unique) > 0) {
        cat("      ", paste(head(pdff_unique, 20), collapse = ", "), "\n")
      }
    }
  }
}

# ==============================================================================
# 7. Figures
# ==============================================================================
cat("\n--- Step 7: Generating figures ---\n")

FIG_DIR <- file.path(BASE_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

if (nrow(coloc_res) > 0) {
  pdf(file.path(FIG_DIR, "broadaway_pdff_coloc_pp4_distribution.pdf"), width = 8, height = 5)
  p1 <- ggplot(coloc_res, aes(x = PP.H4)) +
    geom_histogram(bins = 50, fill = "#2C7BB6", alpha = 0.8) +
    geom_vline(xintercept = c(0.5, 0.8), linetype = "dashed", color = c("orange", "red")) +
    labs(title = "COLOC PP.H4: Broadaway Liver eQTL (N=1,183) x PDFF (N=33,588)",
         subtitle = paste0(n_tested, " genes tested; ",
                          nrow(coloc_res[PP.H4 > 0.5]), " with PP.H4 > 0.5; hg19 native merge"),
         x = "PP.H4 (Posterior Probability of Shared Causal Variant)",
         y = "Count") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 13))
  print(p1)
  dev.off()
  cat("  Saved: broadaway_pdff_coloc_pp4_distribution.pdf\n")

  top30 <- head(coloc_res[PP.H4 > 0.1], 30)
  if (nrow(top30) > 0) {
    top30[, gene := factor(gene, levels = rev(gene))]
    top30[, sig := ifelse(PP.H4 > 0.5, "PP.H4 > 0.5", "PP.H4 <= 0.5")]

    pdf(file.path(FIG_DIR, "broadaway_pdff_coloc_top_genes.pdf"), width = 10, height = 8)
    p2 <- ggplot(top30, aes(x = PP.H4, y = gene, fill = sig)) +
      geom_col(width = 0.7) +
      geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey40") +
      scale_fill_manual(values = c("PP.H4 > 0.5" = "#D7191C", "PP.H4 <= 0.5" = "#2C7BB6")) +
      labs(title = "Top COLOC Hits: Broadaway Liver eQTL x PDFF GWAS (Pazoki 2022)",
           subtitle = "N(eQTL) = 1,183; N(GWAS) = 33,588; hg19 native merge; priors: p12 = 5e-6",
           x = "PP.H4", y = NULL) +
      theme_bw(base_size = 11) +
      theme(plot.title = element_text(face = "bold", size = 13),
            plot.subtitle = element_text(size = 9, color = "grey40"),
            legend.position = "bottom")
    print(p2)
    dev.off()
    cat("  Saved: broadaway_pdff_coloc_top_genes.pdf\n")
  }
}

cat("\n=== Script 35h: Complete ===\n")
cat("End time:", format(Sys.time()), "\n")
