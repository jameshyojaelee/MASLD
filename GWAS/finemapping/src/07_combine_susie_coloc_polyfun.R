#!/usr/bin/env Rscript
# 07_combine_susie_coloc_polyfun.R
# PolyFun UKBB EUR LD-panel analog of 07_combine_susie_coloc{,_topld}.R
# Combines per-chromosome ABF + SuSiE COLOC results across all 17 EUR GWAS
# from the PolyFun re-run at GWAS/finemapping/results/susie_coloc_polyfun/.
#
# Output: GWAS/finemapping/results/susie_coloc_polyfun/gene_level_coloc_polyfun.csv
# Schema: identical to v1 / 1KG / TOP-LD outputs for direct 4-way comparability.
#
# RESEARCH-ONLY (LD-panel sensitivity arm) — NOT the atlas source. This variant
#   processes EUR GWAS only, so it intentionally OMITS the cross-ancestry
#   provenance/tier columns added to the canonical 07 (coloc_*_ancestry,
#   _conf_tier, _ancestry_matched, _headline_cross_anc, _pp4_EUR, n_anc_*): every
#   row here is EUR × EUR-eQTL (ancestry_matched), so those labels are moot. The
#   canonical atlas reads results/susie_coloc/gene_level_coloc.csv (27a:273); do
#   NOT wire this output into 27a/75/217 or any headline claim.
#
# Auto-traverses subdirs of susie_coloc_polyfun/ (no hardcoded GWAS list) so it
# automatically picks up new GWAS results as 8d completes per-chr tasks.
#
# Usage:
#   micromamba run -n rnaseq Rscript GWAS/finemapping/src/07_combine_susie_coloc_polyfun.R

library(data.table)
library(dplyr)

MIN_SNPS <- 100L

BASE_DIR <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FM_DIR <- file.path(BASE_DIR, "GWAS/finemapping")
COLOC_DIR <- file.path(FM_DIR, "results/susie_coloc_polyfun")
LD_FILE_REF <- file.path(FM_DIR, "results/susie_coloc/ld_contamination_clusters.csv")

cat("============================================================\n")
cat("Combining SuSiE+ABF COLOC results (PolyFun UKBB EUR LD)\n")
cat("Input dir :", COLOC_DIR, "\n")
cat("============================================================\n")

if (!dir.exists(COLOC_DIR)) {
  stop("COLOC_DIR does not exist: ", COLOC_DIR, " — has 8d started writing yet?")
}

gwas_dirs <- list.dirs(COLOC_DIR, recursive = FALSE)
all_results <- list()

for (gdir in gwas_dirs) {
  gwas_name <- basename(gdir)
  chr_files <- list.files(gdir, pattern = "susie_coloc_chr.*\\.csv$", full.names = TRUE)
  if (length(chr_files) == 0) {
    chr_files <- list.files(gdir, pattern = "checkpoint_chr.*\\.csv$", full.names = TRUE)
  }
  if (length(chr_files) == 0) {
    cat("  WARNING: No results for", gwas_name, "\n")
    next
  }

  cat("\n--- GWAS:", gwas_name, "---\n")
  gwas_results <- rbindlist(lapply(chr_files, fread), fill = TRUE)
  cat("  Chromosomes:", length(chr_files), "\n")
  cat("  Genes tested:", nrow(gwas_results), "\n")
  cat("  PP.H4.abf > 0.8:", sum(gwas_results$PP.H4.abf > 0.8, na.rm = TRUE), "\n")
  cat("  PP.H4.abf > 0.5:", sum(gwas_results$PP.H4.abf > 0.5, na.rm = TRUE), "\n")
  if ("PP.H4.susie" %in% names(gwas_results)) {
    cat("  PP.H4.susie > 0.8:", sum(gwas_results$PP.H4.susie > 0.8, na.rm = TRUE), "\n")
    cat("  PP.H4.susie > 0.5:", sum(gwas_results$PP.H4.susie > 0.5, na.rm = TRUE), "\n")
  }

  fwrite(gwas_results, file.path(gdir, "susie_coloc_combined.csv"))
  all_results[[length(all_results) + 1]] <- gwas_results
}

if (length(all_results) > 0) {
  master <- rbindlist(all_results, fill = TRUE)

  n_before <- nrow(master)
  master <- master[is.na(n_snps) | n_snps >= MIN_SNPS]
  cat(sprintf("\n  Filtered %d gene-GWAS pairs with n_snps < %d\n", n_before - nrow(master), MIN_SNPS))

  master_file <- file.path(COLOC_DIR, "susie_coloc_all_gwas_polyfun.csv")
  fwrite(master, master_file)
  cat("\n  Master table:", nrow(master), "gene-GWAS entries → ", master_file, "\n")

  # GWAS group assignments (identical mapping; PolyFun is EUR-only so non-EUR
  # rows are absent but the table is kept comprehensive for forward-compat)
  gwas_groups <- data.table(
    gwas_name = c(
      "UKBB_ALT", "UKBB_AST", "UKBB_GGT",
      "2019_31311600_NAFLD_EUR", "2020_32298765_NAFLD_EUR",
      "2021_34841290_NAFLD_EUR",
      "2023_36280732_NAFLD_deCode_EUR", "2023_36280732_NAFLD_Intermountain_EUR",
      "2023_36280732_NAFLD_UKBB_EUR",
      "2021_34128465_PDFF_EUR", "2021_34957434_PDFF_EUR", "2022_36402844_PDFF_EUR",
      "FinnGen_NAFLD", "FinnGen_NASH"
      # cirrhosis/HCC dropped 2026-06-06 (not MASLD-specific): FinnGen_HCC, Ghouse_Cirrhosis, Ghouse_HCC
    ),
    gwas_group = c(
      "UKBB_ALT", "UKBB_AST", "UKBB_GGT",
      "Ghodsian_meta", "Anstee_biopsy",
      "Ghodsian_meta",
      "Sveinbjornsson_deCode", "Sveinbjornsson_Intermountain",
      "Ghodsian_meta",
      "UKBB_PDFF", "UKBB_PDFF", "UKBB_PDFF",
      "FinnGen_R12", "FinnGen_R12"
      # cirrhosis/HCC dropped 2026-06-06 (not MASLD-specific)
    )
  )

  master <- merge(master, gwas_groups, by = "gwas_name", all.x = TRUE)
  master[is.na(gwas_group), gwas_group := gwas_name]

  HAS_SUSIE <- "PP.H4.susie" %in% names(master)

  gene_best <- master %>%
    group_by(gene, ensembl) %>%
    summarise(
      coloc_best_pp4 = max(PP.H4.abf, na.rm = TRUE),
      coloc_best_gwas = if (all(is.na(PP.H4.abf))) NA_character_ else gwas_name[which.max(PP.H4.abf)],
      coloc_n_gwas_h4_05 = sum(PP.H4.abf > 0.5, na.rm = TRUE),
      coloc_n_gwas_h4_08 = sum(PP.H4.abf > 0.8, na.rm = TRUE),
      coloc_n_gwas_tested = n(),
      coloc_n_groups_h4_05 = n_distinct(gwas_group[!is.na(PP.H4.abf) & PP.H4.abf > 0.5]),
      coloc_n_groups_h4_08 = n_distinct(gwas_group[!is.na(PP.H4.abf) & PP.H4.abf > 0.8]),
      .groups = "drop"
    ) %>%
    mutate(coloc_best_pp4 = ifelse(is.infinite(coloc_best_pp4), NA, coloc_best_pp4))

  if (HAS_SUSIE) {
    master_susie <- as.data.frame(master)
    gene_susie <- master_susie %>%
      filter(!is.na(PP.H4.susie)) %>%
      group_by(gene, ensembl) %>%
      summarise(
        coloc_best_susie_pp4       = max(PP.H4.susie, na.rm = TRUE),
        coloc_best_susie_gwas      = gwas_name[which.max(PP.H4.susie)],
        coloc_n_gwas_susie_h4_05   = sum(PP.H4.susie > 0.5, na.rm = TRUE),
        coloc_n_gwas_susie_h4_08   = sum(PP.H4.susie > 0.8, na.rm = TRUE),
        coloc_n_gwas_susie_h4_09   = sum(PP.H4.susie > 0.9, na.rm = TRUE),
        coloc_susie_n_pairs_total  = sum(n_cs_pairs, na.rm = TRUE),
        coloc_susie_success_rate   = mean(method %in% c("susie", "abf_fallback"), na.rm = TRUE),
        coloc_n_groups_susie_h4_05 = n_distinct(gwas_group[!is.na(PP.H4.susie) & PP.H4.susie > 0.5]),
        coloc_n_groups_susie_h4_08 = n_distinct(gwas_group[!is.na(PP.H4.susie) & PP.H4.susie > 0.8]),
        .groups = "drop"
      ) %>%
      mutate(coloc_best_susie_pp4 = ifelse(is.infinite(coloc_best_susie_pp4), NA, coloc_best_susie_pp4))
    gene_best <- left_join(gene_best, gene_susie, by = c("gene", "ensembl"))
    cat(sprintf(
      "  SuSiE PP.H4>0.5: %d | >0.8: %d | >0.9: %d\n",
      sum(gene_best$coloc_n_gwas_susie_h4_05 > 0, na.rm = TRUE),
      sum(gene_best$coloc_n_gwas_susie_h4_08 > 0, na.rm = TRUE),
      sum(gene_best$coloc_n_gwas_susie_h4_09 > 0, na.rm = TRUE)
    ))
  } else {
    for (col in c("coloc_best_susie_pp4", "coloc_best_susie_gwas",
                  "coloc_n_gwas_susie_h4_05", "coloc_n_gwas_susie_h4_08",
                  "coloc_n_gwas_susie_h4_09", "coloc_susie_n_pairs_total",
                  "coloc_susie_success_rate", "coloc_n_groups_susie_h4_05",
                  "coloc_n_groups_susie_h4_08")) gene_best[[col]] <- NA
  }

  # MHC flag
  MHC_CHR <- 6L; MHC_START <- 25000000L; MHC_END <- 35000000L
  mhc_rows <- master[chr == MHC_CHR & !is.na(top_snp)]
  if (nrow(mhc_rows) > 0) {
    mhc_rows[, snp_pos := as.numeric(sub(".*:", "", top_snp))]
    mhc_ensembl <- unique(mhc_rows[snp_pos >= MHC_START & snp_pos <= MHC_END, ensembl])
  } else {
    mhc_ensembl <- character(0)
  }
  gene_best$is_mhc <- gene_best$ensembl %in% mhc_ensembl

  # LD contamination cluster annotation (shared, from v1)
  gene_best <- as.data.table(gene_best)
  if (file.exists(LD_FILE_REF)) {
    ld_clusters <- fread(LD_FILE_REF)
    ld_rows <- list()
    for (r in seq_len(nrow(ld_clusters))) {
      cluster_genes <- trimws(strsplit(ld_clusters$genes[r], ",")[[1]])
      for (g in cluster_genes) {
        ld_rows[[length(ld_rows) + 1]] <- data.table(
          gene = g,
          ld_cluster_best_gene = ld_clusters$best_gene[r],
          ld_cluster_flag = ld_clusters$flag_reason[r]
        )
      }
    }
    ld_dt <- rbindlist(ld_rows)
    ld_dt <- ld_dt[!duplicated(gene)]
    gene_best <- merge(gene_best, ld_dt, by = "gene", all.x = TRUE)
  } else {
    gene_best[, ld_cluster_best_gene := NA_character_]
    gene_best[, ld_cluster_flag := NA_character_]
  }

  gene_file <- file.path(COLOC_DIR, "gene_level_coloc_polyfun.csv")
  fwrite(gene_best, gene_file, na = "NA")
  cat("\nGene-level table:", nrow(gene_best), "genes → ", gene_file, "\n")
  cat("============================================================\n")
}
