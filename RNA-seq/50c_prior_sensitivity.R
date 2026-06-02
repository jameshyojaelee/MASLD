#!/usr/bin/env Rscript
# 50c_prior_sensitivity.R
# ---------------------------------------------------------------------------
# COLOC Prior Sensitivity Analysis (M1)
#
# Re-runs COLOC across a range of p12 priors to assess how sensitive
# results are to the default p12=5e-6 assumption. Compares:
#   - Small-N source: PanUKBB AFR ALT (N=6,636) — expect prior-sensitive
#   - Large-N source: UKBB ALT (N=361,194) — expect robust to prior
#
# Tests p12 = {1e-6, 5e-6, 1e-5, 5e-5, 1e-4}
#
# Output: RNA-seq/results/causal_inference/prior_sensitivity/
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(coloc)
  library(data.table)
  library(rtracklayer)
  library(GenomicRanges)
  library(ggplot2)
})

BASE_DIR    <- Sys.getenv("MASLD_PROJECT_ROOT",
                          "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
EQTL_DIR   <- file.path(BASE_DIR, "data/broadaway_eqtl")
RESULTS_DIR <- file.path(BASE_DIR, "RNA-seq/results/causal_inference/prior_sensitivity")
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

# Sources to test
sources <- list(
  PanUKBB_AFR_ALT = list(
    file = file.path(BASE_DIR, "GWAS/MR_Data/PanUKBB/PanUKBB_AFR_ALT_harmonised_hg38.tsv.gz"),
    N = 6636L, type = "quant", label = "PanUKBB AFR ALT (N=6,636)"
  ),
  UKBB_ALT = list(
    file = file.path(BASE_DIR, "GWAS/MR_Data/GCST90019492_UKBB_ALT_harmonised.tsv.gz"),
    N = 361194L, type = "quant", label = "UKBB ALT (N=361,194)"
  )
)

P12_VALUES <- c(1e-6, 5e-6, 1e-5, 5e-5, 1e-4)
EQTL_N     <- 1183L
COLOC_P1   <- 1e-4
COLOC_P2   <- 1e-4
MIN_MAF    <- 0.01

cat("=== Script 50c: COLOC Prior Sensitivity ===\n")
cat("p12 values:", paste(P12_VALUES, collapse = ", "), "\n")
cat("Start time:", format(Sys.time()), "\n\n")

# ==============================================================================
# 1. Setup
# ==============================================================================
CHAIN_FILE <- file.path(EQTL_DIR, "hg19ToHg38.over.chain")
if (!file.exists(CHAIN_FILE)) {
  CHAIN_GZ <- paste0(CHAIN_FILE, ".gz")
  if (!file.exists(CHAIN_GZ)) {
    download.file(
      "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
      CHAIN_GZ, mode = "wb", quiet = FALSE)
  }
  system2("gunzip", args = c("-k", CHAIN_GZ))
}
chain <- import.chain(CHAIN_FILE)

leads <- fread(file.path(EQTL_DIR, "Liver_eQTL_Meta_Leads_ST3_20240530.tsv"))
all_eGenes <- unique(leads$Gene)
leads[, eqtl_chr := as.integer(sub("_.*", "", Variant))]
gene_chr_map <- unique(leads[!is.na(eqtl_chr), .(Gene, chr = eqtl_chr)])
gene_chr_map <- gene_chr_map[, .N, by = .(Gene, chr)][order(-N)][!duplicated(Gene)]
leads_lookup <- unique(leads[, .(Gene, Ensembl)])

liftover_positions <- function(chr_num, positions) {
  gr <- GRanges(seqnames = paste0("chr", chr_num), ranges = IRanges(start = positions, width = 1))
  lifted <- liftOver(gr, chain)
  n_mapped <- lengths(lifted)
  hg38_pos <- rep(NA_integer_, length(positions))
  idx <- which(n_mapped == 1L)
  if (length(idx) > 0) hg38_pos[idx] <- start(unlist(lifted[idx]))
  hg38_pos
}

chr_eqtl_cache <- list()
load_chr_eqtl <- function(chr_num) {
  key <- as.character(chr_num)
  if (!is.null(chr_eqtl_cache[[key]])) return(chr_eqtl_cache[[key]])
  fname <- file.path(EQTL_DIR, paste0("chr", chr_num, "_marginal_summary_results.tsv"))
  if (!file.exists(fname)) return(NULL)
  dt <- fread(fname)
  dt[, pos_hg38 := liftover_positions(chr_num, POS)]
  dt <- dt[!is.na(pos_hg38)]
  dt[, merge_key := paste0(CHR, ":", pos_hg38)]
  chr_eqtl_cache[[key]] <<- dt
  return(dt)
}

harmonize_alleles <- function(dt) {
  dt <- copy(dt)
  comp <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C")
  dt[, match_type := "none"]
  dt[EA == effect_allele & NEA == other_allele, match_type := "direct"]
  dt[EA == other_allele & NEA == effect_allele, match_type := "flipped"]
  dt[, is_ambiguous := (EA %in% c("A","T") & NEA %in% c("A","T")) |
                        (EA %in% c("C","G") & NEA %in% c("C","G"))]
  dt <- dt[!(is_ambiguous & match_type == "none")]
  dt[match_type == "none" & nchar(EA) == 1 & nchar(NEA) == 1,
     match_type := fifelse(
       comp[EA] == effect_allele & comp[NEA] == other_allele, "direct",
       fifelse(comp[EA] == other_allele & comp[NEA] == effect_allele, "flipped", "none"))]
  dt <- dt[match_type != "none"]
  if (any(dt$match_type == "flipped")) dt[match_type == "flipped", Beta := -Beta]
  dt[, is_ambiguous := NULL]
  return(dt)
}

# ==============================================================================
# 2. Run per source
# ==============================================================================
all_results <- list()

for (src_name in names(sources)) {
  src <- sources[[src_name]]
  cat("\n=== Source:", src$label, "===\n")

  if (!file.exists(src$file)) {
    cat("  File not found, skipping.\n")
    next
  }

  # Load GWAS
  gwas <- fread(src$file)
  gwas[, chr := as.integer(chromosome)]
  gwas[, pos_hg38 := as.integer(base_pair_location)]
  gwas <- gwas[!is.na(chr) & chr %in% 1:22]
  gwas <- gwas[!is.na(beta) & !is.na(standard_error) & standard_error > 0]
  if ("EAF" %in% names(gwas)) {
    gwas[, maf := pmin(EAF, 1 - EAF)]
    gwas <- gwas[is.na(maf) | maf >= MIN_MAF]
  }
  gwas <- gwas[order(chr, pos_hg38, p_value)]
  gwas <- gwas[!duplicated(paste(chr, pos_hg38))]
  gwas[, merge_key := paste0(chr, ":", pos_hg38)]
  cat("  GWAS variants:", format(nrow(gwas), big.mark = ","), "\n")

  # First pass: identify candidate genes from existing results (fast lookup)
  # instead of re-running COLOC for all ~6500 genes
  existing_file <- file.path(BASE_DIR, "RNA-seq/results/causal_inference",
                              tolower(gsub(" .*", "", gsub("PanUKBB_", "panukbb_", src_name))),
                              "coloc_results.csv")
  # Map source names to result directories
  result_dir_map <- c(
    PanUKBB_AFR_ALT = "panukbb_afr_alt",
    UKBB_ALT = ""  # handled separately below
  )
  if (src_name == "PanUKBB_AFR_ALT") {
    existing_file <- file.path(BASE_DIR, "RNA-seq/results/causal_inference/panukbb_afr_alt/coloc_results.csv")
  } else if (src_name == "UKBB_ALT") {
    # UKBB ALT uses broadaway COLOC results (ALT = broadaway_coloc_pp4)
    existing_file <- file.path(BASE_DIR, "RNA-seq/results/causal_inference/broadaway_ukbb/coloc_results.csv")
  }

  cat("  Identifying candidate genes (PP.H4 > 0.3) from existing results...\n")
  candidate_genes <- character()
  if (file.exists(existing_file)) {
    existing_dt <- fread(existing_file)
    candidate_genes <- existing_dt[PP.H4 > 0.3, unique(gene)]
    cat("  Found", length(candidate_genes), "candidates from", basename(existing_file), "\n")
  } else {
    cat("  No existing results at", existing_file, "— running first pass from scratch\n")
    for (gene_name in all_eGenes) {
      chr_info <- gene_chr_map[Gene == gene_name]
      if (nrow(chr_info) == 0) next

      eqtl_data <- load_chr_eqtl(chr_info$chr[1])
      if (is.null(eqtl_data)) next

      gene_eqtl <- eqtl_data[GeneSymbol == gene_name]
      if (nrow(gene_eqtl) == 0) {
        gene_ensg <- leads_lookup$Ensembl[leads_lookup$Gene == gene_name]
        if (length(gene_ensg) > 0) gene_eqtl <- eqtl_data[ENSG %in% gene_ensg]
      }
      if (nrow(gene_eqtl) < 10) next

      merged <- merge(gene_eqtl, gwas, by = "merge_key", suffixes = c(".eqtl", ".gwas"))
      merged <- merged[!is.na(Beta) & !is.na(SE) & SE > 0 &
                       !is.na(beta) & !is.na(standard_error) & standard_error > 0]
      if (nrow(merged) < 10) next

      merged <- harmonize_alleles(merged)
      if (nrow(merged) < 10) next

      d1 <- list(snp = merged$merge_key, beta = merged$beta,
                 varbeta = (merged$standard_error)^2,
                 type = "quant", sdY = 1, N = src$N)
      d2 <- list(snp = merged$merge_key, beta = merged$Beta,
                 varbeta = (merged$SE)^2,
                 type = "quant", sdY = 1, N = EQTL_N)

      res <- tryCatch(coloc.abf(d1, d2, p1 = COLOC_P1, p2 = COLOC_P2, p12 = 5e-6),
                      error = function(e) NULL)

      if (!is.null(res) && res$summary["PP.H4.abf"] > 0.3) {
        candidate_genes <- c(candidate_genes, gene_name)
      }
    }
  }
  cat("  Candidates with PP.H4 > 0.3:", length(candidate_genes), "\n")

  # Second pass: re-run candidates across all p12 values
  cat("  Running sensitivity across p12 values...\n")

  for (gene_name in candidate_genes) {
    chr_info <- gene_chr_map[Gene == gene_name]
    eqtl_data <- load_chr_eqtl(chr_info$chr[1])
    gene_eqtl <- eqtl_data[GeneSymbol == gene_name]
    if (nrow(gene_eqtl) == 0) {
      gene_ensg <- leads_lookup$Ensembl[leads_lookup$Gene == gene_name]
      if (length(gene_ensg) > 0) gene_eqtl <- eqtl_data[ENSG %in% gene_ensg]
    }

    merged <- merge(gene_eqtl, gwas, by = "merge_key", suffixes = c(".eqtl", ".gwas"))
    merged <- merged[!is.na(Beta) & !is.na(SE) & SE > 0 &
                     !is.na(beta) & !is.na(standard_error) & standard_error > 0]
    merged <- harmonize_alleles(merged)
    if (nrow(merged) < 10) next

    d1 <- list(snp = merged$merge_key, beta = merged$beta,
               varbeta = (merged$standard_error)^2,
               type = "quant", sdY = 1, N = src$N)
    d2 <- list(snp = merged$merge_key, beta = merged$Beta,
               varbeta = (merged$SE)^2,
               type = "quant", sdY = 1, N = EQTL_N)

    for (p12 in P12_VALUES) {
      res <- tryCatch(coloc.abf(d1, d2, p1 = COLOC_P1, p2 = COLOC_P2, p12 = p12),
                      error = function(e) NULL)
      if (!is.null(res)) {
        all_results[[length(all_results) + 1]] <- data.table(
          source = src_name, source_N = src$N,
          gene = gene_name, p12 = p12,
          PP.H4 = res$summary["PP.H4.abf"],
          PP.H3 = res$summary["PP.H3.abf"],
          n_snps = nrow(merged)
        )
      }
    }
  }
}

# ==============================================================================
# 3. Analyze + plot
# ==============================================================================
cat("\n--- Analysis ---\n")

if (length(all_results) == 0) {
  cat("  No results.\n")
  quit(save = "no", status = 0)
}

sens <- rbindlist(all_results)
fwrite(sens, file.path(RESULTS_DIR, "prior_sensitivity_results.csv"))

# Per-source summary: n genes PP.H4 > 0.5 at each p12
cat("\n  Genes with PP.H4 > 0.5 by p12:\n")
summary_tab <- sens[, .(n_sig = sum(PP.H4 > 0.5)), by = .(source, p12)]
summary_tab <- dcast(summary_tab, source ~ p12, value.var = "n_sig")
print(summary_tab)
fwrite(summary_tab, file.path(RESULTS_DIR, "prior_sensitivity_summary.csv"))

# Fragile genes: PP.H4 crosses 0.5 threshold across p12 range
fragile <- sens[, .(min_pp4 = min(PP.H4), max_pp4 = max(PP.H4),
                     pp4_at_default = PP.H4[p12 == 5e-6],
                     pp4_range = max(PP.H4) - min(PP.H4)),
                by = .(source, gene)]
fragile[, crosses_threshold := (min_pp4 < 0.5) & (max_pp4 > 0.5)]
fragile_genes <- fragile[crosses_threshold == TRUE]
setorder(fragile_genes, source, -pp4_range)
fwrite(fragile_genes, file.path(RESULTS_DIR, "prior_fragile_genes.csv"))

cat("\n  Fragile genes (cross PP.H4=0.5 threshold across p12 range):\n")
for (src_name in unique(fragile_genes$source)) {
  fg <- fragile_genes[source == src_name]
  cat("    ", src_name, ":", nrow(fg), "fragile genes\n")
  for (i in seq_len(min(5, nrow(fg)))) {
    cat("      ", fg$gene[i], ": range=", round(fg$min_pp4[i], 3), "-",
        round(fg$max_pp4[i], 3), "\n")
  }
}

# Figure: sensitivity curves
pdf(file.path(RESULTS_DIR, "prior_sensitivity_curves.pdf"), width = 10, height = 6)
p <- ggplot(sens[, .(n_sig = sum(PP.H4 > 0.5)), by = .(source, p12)],
            aes(x = p12, y = n_sig, color = source)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_vline(xintercept = 5e-6, linetype = "dashed", color = "grey50") +
  scale_x_log10(labels = scales::scientific) +
  labs(title = "COLOC Prior Sensitivity: p12 vs Number of Hits",
       subtitle = "Dashed line = default p12=5e-6",
       x = "p12 prior", y = "Genes with PP.H4 > 0.5",
       color = "GWAS Source") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom")
print(p)
dev.off()
cat("  Saved: prior_sensitivity_curves.pdf\n")

# Per-gene PP.H4 heatmap for fragile genes
if (nrow(fragile_genes) > 0) {
  fragile_data <- sens[gene %in% fragile_genes$gene]
  fragile_data[, gene := factor(gene, levels = rev(unique(fragile_genes$gene)))]
  fragile_data[, p12_label := sprintf("%.0e", p12)]

  pdf(file.path(RESULTS_DIR, "prior_sensitivity_fragile_heatmap.pdf"),
      width = 10, height = max(4, nrow(fragile_genes) * 0.3 + 2))
  p2 <- ggplot(fragile_data, aes(x = p12_label, y = gene, fill = PP.H4)) +
    geom_tile() +
    geom_text(aes(label = round(PP.H4, 2)), size = 2.5) +
    scale_fill_gradient2(low = "white", mid = "#FEE08B", high = "#D73027",
                         midpoint = 0.5) +
    facet_wrap(~source, scales = "free_y") +
    labs(title = "PP.H4 Across p12 Priors: Fragile Genes",
         x = "p12 prior", y = NULL) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold"))
  print(p2)
  dev.off()
  cat("  Saved: prior_sensitivity_fragile_heatmap.pdf\n")
}

cat("\n=== Script 50c: Complete ===\n")
cat("End time:", format(Sys.time()), "\n")
