#!/usr/bin/env Rscript
# 35x_gtex_coloc_ukbb.R
# ---------------------------------------------------------------------------
# Colocalization: GTEx v8 Liver eQTL (N=208) x UKBB ALT GWAS (N=343,850)
#
# Independent replication of the Broadaway COLOC (Script 35b) using the
# smaller but widely-used GTEx v8 Liver eQTL dataset. Both GWAS and eQTL
# are on hg38, so no liftover is needed.
#
# Key differences from Script 35b (Broadaway):
#   - eQTL source: GTEx v8 Liver (N=208) vs Broadaway (N=1,183)
#   - No liftover needed (both GTEx and UKBB are hg38)
#   - Gene IDs are Ensembl (ENSG...) — mapped to symbols via GENCODE cache
#   - Larger file (163M rows) — processed chromosome-by-chromosome
#
# Inputs:
#   - GWAS/MR_Data/GTEx_v8_Liver_eQTL.tsv.gz (3.1 GB, 163M rows)
#   - GWAS/MR_Data/GCST90019492_UKBB_ALT_harmonised.tsv.gz
#   - gene annotation cache: human_ensg_to_symbol.tsv
#
# Outputs:
#   - results/causal_inference/gtex_v8_ukbb_alt/coloc_results.csv
#   - results/causal_inference/gtex_v8_ukbb_alt/coloc_summary.csv
#   - results/causal_inference/gtex_v8_ukbb_alt/coloc_sensitivity/ (PDFs)
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(coloc)
  library(data.table)
  library(ggplot2)
})

cat("=== Script 35x: GTEx v8 Liver eQTL COLOC x UKBB ALT GWAS ===\n")
cat("Start time:", format(Sys.time()), "\n\n")

# ==============================================================================
# Configuration
# ==============================================================================
BASE_DIR    <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
EQTL_FILE   <- file.path(BASE_DIR, "GWAS/MR_Data/GTEx_v8_Liver_eQTL.tsv.gz")
GWAS_FILE   <- file.path(BASE_DIR, "GWAS/MR_Data/GCST90019492_UKBB_ALT_harmonised.tsv.gz")
GENE_CACHE  <- file.path(BASE_DIR,
  "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/gene_annotation/human_ensg_to_symbol.tsv")
RESULTS_DIR <- file.path(BASE_DIR, "RNA-seq/results/causal_inference/gtex_v8_ukbb_alt")
SENS_DIR    <- file.path(RESULTS_DIR, "coloc_sensitivity")

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SENS_DIR, recursive = TRUE, showWarnings = FALSE)

# GWAS parameters (UKBB ALT: quantitative, N=343,850)
GWAS_N    <- 343850L
GWAS_TYPE <- "quant"

# GTEx v8 Liver eQTL parameters
EQTL_N    <- 208L

# COLOC priors (same as 35b)
COLOC_P1  <- 1e-4
COLOC_P2  <- 1e-4
COLOC_P12 <- 5e-6

# Minimum shared SNPs for COLOC
MIN_SNPS  <- 10L

# ==============================================================================
# 1. Load gene symbol mapping
# ==============================================================================
cat("--- Step 1: Loading gene annotation cache ---\n")

if (file.exists(GENE_CACHE)) {
  gene_annot <- fread(GENE_CACHE)
  gene_annot <- gene_annot[!duplicated(gene_base), .(ensembl_id = gene_base, symbol)]
  cat("  Gene annotations loaded:", nrow(gene_annot), "mappings\n")
} else {
  cat("  WARNING: Gene annotation cache not found at", GENE_CACHE, "\n")
  cat("  Will report Ensembl IDs only (no gene symbols)\n")
  gene_annot <- data.table(ensembl_id = character(), symbol = character())
}

# Helper: map Ensembl ID to symbol
ensembl_to_symbol <- function(ensg_id) {
  sym <- gene_annot$symbol[match(ensg_id, gene_annot$ensembl_id)]
  ifelse(is.na(sym) | sym == "", ensg_id, sym)
}

# ==============================================================================
# 2. Load UKBB ALT GWAS
# ==============================================================================
cat("\n--- Step 2: Loading UKBB ALT GWAS (GCST90019492) ---\n")

gwas <- fread(GWAS_FILE)
cat("  Raw GWAS rows:", nrow(gwas), "\n")

# Filter to valid autosomal rows with beta+SE
gwas <- gwas[!is.na(chromosome) & !is.na(base_pair_location) & !is.na(p_value)]
gwas[, chr := as.integer(chromosome)]
gwas[, pos := as.integer(base_pair_location)]
gwas <- gwas[!is.na(chr) & chr %in% 1:22]
gwas <- gwas[!is.na(beta) & !is.na(standard_error) & standard_error > 0]

# Deduplicate by position (keep lowest p-value)
gwas <- gwas[order(chr, pos, p_value)]
gwas <- gwas[!duplicated(paste(chr, pos))]
gwas[, merge_key := paste0(chr, ":", pos)]

cat("  GWAS after dedup:", nrow(gwas), "variants\n")
cat("  Genome-wide significant (p<5e-8):", sum(gwas$p_value < 5e-8), "\n")

# Index GWAS by chromosome for fast per-chr lookups
setkey(gwas, chr, pos)

# ==============================================================================
# 3. Load GTEx v8 Liver eQTL (chromosome by chromosome)
# ==============================================================================
cat("\n--- Step 3: Loading GTEx v8 Liver eQTL ---\n")
cat("  File:", EQTL_FILE, "\n")
cat("  Reading with column selection for memory efficiency...\n")

# Load full file with selected columns
# GTEx columns: variant, r2, pvalue, molecular_trait_object_id, molecular_trait_id,
#               maf, gene_id, median_tpm, beta, se, an, ac, chromosome, position,
#               ref, alt, type, rsid
eqtl_all <- fread(
  EQTL_FILE,
  select = c("chromosome", "position", "ref", "alt", "beta", "se",
             "pvalue", "maf", "gene_id", "rsid")
)
cat("  Raw eQTL rows:", nrow(eqtl_all), "\n")

# Filter to valid data
eqtl_all <- eqtl_all[!is.na(chromosome) & !is.na(position) &
                      !is.na(beta) & !is.na(se) & se > 0]
eqtl_all[, chr := as.integer(chromosome)]
eqtl_all <- eqtl_all[!is.na(chr) & chr %in% 1:22]

# Create merge key (hg38 — same build as UKBB)
eqtl_all[, merge_key := paste0(chr, ":", position)]

n_genes <- uniqueN(eqtl_all$gene_id)
cat("  eQTL after filtering:", nrow(eqtl_all), "rows\n")
cat("  Unique eGenes:", n_genes, "\n")
cat("  Chromosomes:", paste(sort(unique(eqtl_all$chr)), collapse = ", "), "\n")

# ==============================================================================
# 4. Helper functions
# ==============================================================================

# Allele harmonization (GTEx ref/alt vs UKBB effect_allele/other_allele)
harmonize_alleles <- function(dt) {
  dt <- copy(dt)
  comp <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C")

  dt[, match_type := "none"]

  # GTEx: ref (reference), alt (alternate/effect allele for beta)
  # UKBB: effect_allele, other_allele
  # Direct: alt == effect_allele AND ref == other_allele
  dt[alt == effect_allele & ref == other_allele, match_type := "direct"]
  # Flipped: alt == other_allele AND ref == effect_allele
  dt[alt == other_allele & ref == effect_allele, match_type := "flipped"]

  # Remove strand-ambiguous variants that failed direct/flipped match
  dt[, is_ambiguous := (ref %in% c("A","T") & alt %in% c("A","T")) |
                        (ref %in% c("C","G") & alt %in% c("C","G"))]
  n_ambig_unmatched <- sum(dt$is_ambiguous & dt$match_type == "none", na.rm = TRUE)
  dt <- dt[!(is_ambiguous & match_type == "none")]

  # Complement matching for remaining non-ambiguous single-base SNPs
  dt[match_type == "none" & nchar(ref) == 1 & nchar(alt) == 1,
     match_type := fifelse(
       comp[alt] == effect_allele & comp[ref] == other_allele, "direct",
       fifelse(comp[alt] == other_allele & comp[ref] == effect_allele, "flipped",
               "none")
     )]

  n_direct  <- sum(dt$match_type == "direct")
  n_flipped <- sum(dt$match_type == "flipped")
  n_none    <- sum(dt$match_type == "none")

  dt <- dt[match_type != "none"]

  # Flip eQTL beta when alleles are swapped
  if (any(dt$match_type == "flipped")) {
    dt[match_type == "flipped", beta.eqtl := -beta.eqtl]
  }

  dt[, is_ambiguous := NULL]
  return(dt)
}

# ==============================================================================
# 5. Run COLOC per eGene (chromosome by chromosome)
# ==============================================================================
cat("\n--- Step 5: Running COLOC ---\n")

all_genes <- unique(eqtl_all$gene_id)
cat("  Total eGenes to test:", length(all_genes), "\n\n")

results_list <- list()
n_tested     <- 0L
n_skipped    <- 0L
n_low_snps   <- 0L
n_no_harm    <- 0L
n_error      <- 0L

# Process chromosome by chromosome for memory-efficient GWAS subsetting
for (chr_num in sort(unique(eqtl_all$chr))) {
  cat("  Processing chromosome", chr_num, "...\n")

  # Subset eQTL and GWAS for this chromosome
  eqtl_chr <- eqtl_all[chr == chr_num]
  gwas_chr <- gwas[chr == chr_num]

  if (nrow(eqtl_chr) == 0 || nrow(gwas_chr) == 0) {
    cat("    Skipping: no data\n")
    next
  }

  chr_genes <- unique(eqtl_chr$gene_id)
  cat("    eGenes on chr", chr_num, ":", length(chr_genes), "\n")

  for (gid in chr_genes) {
    # Use explicit variable to avoid data.table column/variable name conflict
    gene_eqtl <- eqtl_chr[gene_id == gid]

    if (nrow(gene_eqtl) < MIN_SNPS) {
      n_low_snps <- n_low_snps + 1L
      next
    }

    # Merge eQTL with GWAS on merge_key (hg38, same build)
    merged <- merge(
      gene_eqtl[, .(merge_key, ref, alt, beta.eqtl = beta, se.eqtl = se,
                     pvalue.eqtl = pvalue, maf, rsid, position)],
      gwas_chr[, .(merge_key, effect_allele, other_allele,
                   beta.gwas = beta, se.gwas = standard_error, p_value)],
      by = "merge_key"
    )

    if (nrow(merged) < MIN_SNPS) {
      n_low_snps <- n_low_snps + 1L
      next
    }

    # Filter to rows with valid beta/SE in both datasets
    merged <- merged[!is.na(beta.eqtl) & !is.na(se.eqtl) & se.eqtl > 0 &
                     !is.na(beta.gwas) & !is.na(se.gwas) & se.gwas > 0]

    if (nrow(merged) < MIN_SNPS) {
      n_low_snps <- n_low_snps + 1L
      next
    }

    # Allele harmonization
    merged <- harmonize_alleles(merged)

    if (nrow(merged) < MIN_SNPS) {
      n_no_harm <- n_no_harm + 1L
      next
    }

    # Deduplicate by merge_key (keep lowest eQTL p-value)
    merged <- merged[order(pvalue.eqtl)]
    merged <- merged[!duplicated(merge_key)]

    if (nrow(merged) < MIN_SNPS) {
      n_low_snps <- n_low_snps + 1L
      next
    }

    # Prepare COLOC datasets with beta/varbeta
    dataset1 <- list(
      snp     = merged$merge_key,
      beta    = merged$beta.gwas,
      varbeta = (merged$se.gwas)^2,
      type    = GWAS_TYPE,
      sdY     = 1,
      N       = GWAS_N
    )

    dataset2 <- list(
      snp     = merged$merge_key,
      beta    = merged$beta.eqtl,
      varbeta = (merged$se.eqtl)^2,
      type    = "quant",
      sdY     = 1,
      N       = EQTL_N
    )

    my.res <- tryCatch(
      coloc.abf(dataset1, dataset2, p1 = COLOC_P1, p2 = COLOC_P2, p12 = COLOC_P12),
      error = function(e) {
        n_error <<- n_error + 1L
        NULL
      }
    )

    if (!is.null(my.res)) {
      n_tested <- n_tested + 1L

      # Extract top SNP
      top_snp_row <- my.res$results[which.max(my.res$results$SNP.PP.H4), ]

      # Map Ensembl to symbol
      sym <- ensembl_to_symbol(gid)

      results_list[[gid]] <- data.table(
        gene        = gid,
        symbol      = sym,
        chr         = chr_num,
        PP.H0       = my.res$summary["PP.H0.abf"],
        PP.H1       = my.res$summary["PP.H1.abf"],
        PP.H2       = my.res$summary["PP.H2.abf"],
        PP.H3       = my.res$summary["PP.H3.abf"],
        PP.H4       = my.res$summary["PP.H4.abf"],
        n_snps      = nrow(merged),
        n_eqtl_snps = nrow(gene_eqtl),
        top_snp     = top_snp_row$snp,
        top_snp_PP  = top_snp_row$SNP.PP.H4,
        eqtl_source = "GTEx_v8_Liver",
        eqtl_N      = EQTL_N,
        gwas_source = "UKBB_ALT",
        gwas_N      = GWAS_N,
        method      = "abf_beta_varbeta"
      )

      # Sensitivity analysis for promising hits
      if (my.res$summary["PP.H4.abf"] > 0.3) {
        safe_name <- gsub("[^A-Za-z0-9_]", "_", gid)
        tryCatch({
          pdf(file.path(SENS_DIR, paste0("sensitivity_", safe_name, ".pdf")),
              width = 10, height = 8)
          sensitivity(my.res, rule = "H4 > 0.5")
          dev.off()
        }, error = function(e) {
          tryCatch(dev.off(), error = function(e2) NULL)
        })
      }

      if (n_tested %% 500 == 0) {
        cat("    [", n_tested, "tested,", n_low_snps, "low-SNP,",
            n_no_harm, "no-harmonize,", n_skipped, "skipped,",
            n_error, "errors]\n")
      }
    }
  }

  cat("    chr", chr_num, "done. Cumulative tested:", n_tested, "\n")
}

# ==============================================================================
# 6. Compile and save results
# ==============================================================================
cat("\n--- Step 6: Compiling results ---\n")

cat("  Genes tested:", n_tested, "\n")
cat("  Genes skipped (< 10 shared SNPs):", n_low_snps, "\n")
cat("  Genes skipped (allele harmonization):", n_no_harm, "\n")
cat("  Genes with COLOC errors:", n_error, "\n")

if (length(results_list) == 0) {
  cat("  No COLOC results.\n")
  fwrite(data.table(gene = character(), symbol = character(), PP.H4 = numeric()),
         file.path(RESULTS_DIR, "coloc_results.csv"))
  quit(save = "no", status = 0)
}

coloc_res <- rbindlist(results_list, fill = TRUE)
setorder(coloc_res, -PP.H4)

fwrite(coloc_res, file.path(RESULTS_DIR, "coloc_results.csv"))

cat("\n  COLOC Summary (GTEx v8 Liver N=208 x UKBB ALT N=343,850):\n")
cat("    Genes tested:", n_tested, "\n")
cat("    PP.H4 > 0.8:", nrow(coloc_res[PP.H4 > 0.8]), "\n")
cat("    PP.H4 > 0.5:", nrow(coloc_res[PP.H4 > 0.5]), "\n")
cat("    PP.H4 > 0.3:", nrow(coloc_res[PP.H4 > 0.3]), "\n")

# Top hits
top_hits <- coloc_res[PP.H4 > 0.3]
if (nrow(top_hits) > 0) {
  cat("\n  Top COLOC hits (PP.H4 > 0.3):\n")
  for (i in seq_len(min(nrow(top_hits), 30))) {
    cat("    ", top_hits$symbol[i], " (", top_hits$gene[i], "): PP.H4 =",
        round(top_hits$PP.H4[i], 3),
        "(", top_hits$n_snps[i], "SNPs, top:", top_hits$top_snp[i], ")\n")
  }
}

# Save summary
summary_dt <- data.table(
  metric = c("genes_tested", "genes_skipped_low_snps",
             "genes_skipped_harmonization", "genes_with_errors",
             "PP.H4_gt_0.8", "PP.H4_gt_0.5", "PP.H4_gt_0.3",
             "eqtl_source", "eqtl_N", "gwas_source", "gwas_N"),
  value = c(n_tested, n_low_snps, n_no_harm, n_error,
            nrow(coloc_res[PP.H4 > 0.8]),
            nrow(coloc_res[PP.H4 > 0.5]),
            nrow(coloc_res[PP.H4 > 0.3]),
            "GTEx_v8_Liver", EQTL_N,
            "UKBB_ALT_GCST90019492", GWAS_N)
)
fwrite(summary_dt, file.path(RESULTS_DIR, "coloc_summary.csv"))

# ==============================================================================
# 7. Validation: Compare with Broadaway COLOC (Script 35b)
# ==============================================================================
cat("\n--- Step 7: Validation ---\n")

# Cross-reference with Broadaway results
broadaway_file <- file.path(BASE_DIR, "RNA-seq/results/causal_inference/broadaway_ukbb/coloc_results.csv")
if (file.exists(broadaway_file)) {
  broadaway_res <- fread(broadaway_file)
  if (nrow(broadaway_res) > 0 && "PP.H4" %in% names(broadaway_res)) {
    cat("\n  Comparison with Script 35b (Broadaway N=1,183 x UKBB ALT):\n")
    cat("    Broadaway genes PP.H4 > 0.5:", nrow(broadaway_res[PP.H4 > 0.5]), "\n")
    cat("    GTEx v8 genes PP.H4 > 0.5:", nrow(coloc_res[PP.H4 > 0.5]), "\n")

    # Broadaway uses gene symbols; GTEx uses Ensembl. Match on both.
    # Broadaway has gene (symbol) and ensembl columns
    if ("ensembl" %in% names(broadaway_res)) {
      # Clean Ensembl IDs (strip version, handle semicolons)
      broadaway_res[, ensembl_clean := sub(";.*", "", sub("\\..*", "", ensembl))]
      comp <- merge(
        broadaway_res[, .(ensembl_clean, gene_broadaway = gene,
                          PP.H4_broadaway = PP.H4)],
        coloc_res[, .(ensembl_clean = gene, symbol,
                      PP.H4_gtex = PP.H4)],
        by = "ensembl_clean"
      )
    } else {
      # Fallback: match via symbol
      comp <- merge(
        broadaway_res[, .(gene, PP.H4_broadaway = PP.H4)],
        coloc_res[, .(gene = symbol, PP.H4_gtex = PP.H4)],
        by = "gene"
      )
    }

    cat("    Shared genes tested:", nrow(comp), "\n")

    if (nrow(comp) > 5) {
      rho <- cor(comp$PP.H4_broadaway, comp$PP.H4_gtex,
                 use = "complete.obs", method = "spearman")
      r   <- cor(comp$PP.H4_broadaway, comp$PP.H4_gtex,
                 use = "complete.obs", method = "pearson")
      cat("    Spearman rho (PP.H4):", round(rho, 3), "\n")
      cat("    Pearson r (PP.H4):", round(r, 3), "\n")

      both_sig <- comp[PP.H4_broadaway > 0.5 & PP.H4_gtex > 0.5]
      broadaway_only <- comp[PP.H4_broadaway > 0.5 & PP.H4_gtex <= 0.5]
      gtex_only <- comp[PP.H4_broadaway <= 0.5 & PP.H4_gtex > 0.5]
      cat("    Both PP.H4 > 0.5:", nrow(both_sig), "\n")
      cat("    Broadaway-only PP.H4 > 0.5:", nrow(broadaway_only), "\n")
      cat("    GTEx-only PP.H4 > 0.5:", nrow(gtex_only), "\n")

      if (nrow(both_sig) > 0) {
        cat("    Concordant hits:\n")
        both_sig <- both_sig[order(-PP.H4_gtex)]
        for (i in seq_len(min(nrow(both_sig), 20))) {
          label <- if ("gene_broadaway" %in% names(both_sig)) {
            both_sig$gene_broadaway[i]
          } else {
            both_sig$gene[i]
          }
          cat("      ", label, ": Broadaway=", round(both_sig$PP.H4_broadaway[i], 3),
              " GTEx=", round(both_sig$PP.H4_gtex[i], 3), "\n")
        }
      }
    }
  }
} else {
  cat("  Broadaway COLOC results not found — skipping comparison.\n")
}

# Zenodo NAFLD gene validation (same targets as 35b)
zenodo_genes_sym <- c("HSD17B13", "EFHD1", "CHEK2", "MAP1LC3A", "SPTLC3",
                       "APOH", "MFHAS1", "ETS2", "LCOR", "USP40", "ALAD")
cat("\n  Zenodo NAFLD gene validation (by symbol):\n")
for (g in zenodo_genes_sym) {
  hit <- coloc_res[symbol == g]
  if (nrow(hit) > 0) {
    cat("    ", g, ": PP.H4 =", round(hit$PP.H4[1], 3),
        "(", hit$n_snps[1], "SNPs)\n")
  } else {
    cat("    ", g, ": not tested (no GTEx eQTL or insufficient overlap)\n")
  }
}

# Known drug targets
drug_targets <- c("THRB", "NR1H4", "PPARA", "PPARG", "PPARD")
cat("\n  Known drug target validation (by symbol):\n")
for (g in drug_targets) {
  hit <- coloc_res[symbol == g]
  if (nrow(hit) > 0) {
    cat("    ", g, ": PP.H4 =", round(hit$PP.H4[1], 3),
        "(", hit$n_snps[1], "SNPs)\n")
  } else {
    cat("    ", g, ": not tested\n")
  }
}

# ==============================================================================
# 8. Figures
# ==============================================================================
cat("\n--- Step 8: Generating figures ---\n")

FIG_DIR <- file.path(RESULTS_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

if (nrow(coloc_res) > 0) {
  # PP.H4 distribution
  pdf(file.path(FIG_DIR, "gtex_ukbb_coloc_pp4_distribution.pdf"), width = 8, height = 5)
  p1 <- ggplot(coloc_res, aes(x = PP.H4)) +
    geom_histogram(bins = 50, fill = "#4DAF4A", alpha = 0.8) +
    geom_vline(xintercept = c(0.5, 0.8), linetype = "dashed", color = c("orange", "red")) +
    labs(title = "COLOC PP.H4: GTEx v8 Liver eQTL (N=208) x UKBB ALT (N=343,850)",
         subtitle = paste0(n_tested, " genes tested; ",
                          nrow(coloc_res[PP.H4 > 0.5]), " with PP.H4 > 0.5"),
         x = "PP.H4 (Posterior Probability of Shared Causal Variant)",
         y = "Count") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 13))
  print(p1)
  dev.off()
  cat("  Saved: gtex_ukbb_coloc_pp4_distribution.pdf\n")

  # Top genes barplot
  top30 <- head(coloc_res[PP.H4 > 0.1], 30)
  if (nrow(top30) > 0) {
    top30[, label := ifelse(symbol != gene, symbol, gene)]
    top30[, label := factor(label, levels = rev(label))]
    top30[, sig := ifelse(PP.H4 > 0.5, "PP.H4 > 0.5", "PP.H4 <= 0.5")]

    pdf(file.path(FIG_DIR, "gtex_ukbb_coloc_top_genes.pdf"), width = 10, height = 8)
    p2 <- ggplot(top30, aes(x = PP.H4, y = label, fill = sig)) +
      geom_col(width = 0.7) +
      geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey40") +
      scale_fill_manual(values = c("PP.H4 > 0.5" = "#D7191C", "PP.H4 <= 0.5" = "#4DAF4A")) +
      labs(title = "Top COLOC Hits: GTEx v8 Liver eQTL x UKBB ALT GWAS",
           subtitle = paste0("N(eQTL) = 208; N(GWAS) = 343,850; priors: p12 = 5e-6; beta/varbeta"),
           x = "PP.H4", y = NULL) +
      theme_bw(base_size = 11) +
      theme(plot.title = element_text(face = "bold", size = 13),
            plot.subtitle = element_text(size = 9, color = "grey40"),
            legend.position = "bottom")
    print(p2)
    dev.off()
    cat("  Saved: gtex_ukbb_coloc_top_genes.pdf\n")
  }

  # Scatter: GTEx vs Broadaway PP.H4 (if comparison data exists)
  if (exists("comp") && nrow(comp) > 5) {
    pdf(file.path(FIG_DIR, "gtex_vs_broadaway_pp4_scatter.pdf"), width = 7, height = 7)
    p3 <- ggplot(comp, aes(x = PP.H4_broadaway, y = PP.H4_gtex)) +
      geom_point(alpha = 0.4, size = 1.5) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
      geom_hline(yintercept = 0.5, linetype = "dotted", color = "grey50") +
      geom_vline(xintercept = 0.5, linetype = "dotted", color = "grey50") +
      labs(title = "COLOC PP.H4 Comparison: Broadaway vs GTEx v8 Liver",
           subtitle = paste0(nrow(comp), " shared genes; Spearman rho = ",
                            round(cor(comp$PP.H4_broadaway, comp$PP.H4_gtex,
                                      use = "complete.obs", method = "spearman"), 3)),
           x = "PP.H4 (Broadaway, N=1,183)",
           y = "PP.H4 (GTEx v8, N=208)") +
      theme_bw(base_size = 11) +
      theme(plot.title = element_text(face = "bold", size = 13))
    print(p3)
    dev.off()
    cat("  Saved: gtex_vs_broadaway_pp4_scatter.pdf\n")
  }
}

cat("\n=== Script 35x: Complete ===\n")
cat("End time:", format(Sys.time()), "\n")
cat("Results saved to:", RESULTS_DIR, "\n")
