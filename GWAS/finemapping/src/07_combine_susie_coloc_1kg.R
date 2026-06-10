#!/usr/bin/env Rscript
# 07_combine_susie_coloc_1kg.R
# 1KG EUR LD-panel analog of 07_combine_susie_coloc.R
# Combines per-chromosome ABF + SuSiE COLOC results across all 17 EUR GWAS
# from the 1KG-LD re-run at GWAS/finemapping/results/susie_coloc_1kg/.
#
# Output: GWAS/finemapping/results/susie_coloc_1kg/gene_level_coloc_1kg.csv
# Schema: identical to v1 gene_level_coloc.csv for direct comparability.
#
# Differences vs 07_combine_susie_coloc.R:
#   - COLOC_DIR points at susie_coloc_1kg/
#   - Per-GWAS combined CSV name unchanged (susie_coloc_combined.csv)
#   - Master / gene-level filenames suffixed _1kg
#   - LD contamination cluster annotation is read from the v1 dir
#     (susie_coloc/ld_contamination_clusters.csv) — same gene-level
#     annotation applies regardless of LD panel; this keeps the
#     ld_cluster_* columns identical between v1 and 1KG outputs and
#     does NOT write into susie_coloc/.
#
# Usage:
#   micromamba run -n rnaseq Rscript GWAS/finemapping/src/07_combine_susie_coloc_1kg.R

library(data.table)
library(dplyr)

MIN_SNPS <- 100L  # Post-hoc filter: exclude gene-GWAS pairs with < 100 overlapping SNPs

BASE_DIR <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FM_DIR <- file.path(BASE_DIR, "GWAS/finemapping")
COLOC_DIR <- file.path(FM_DIR, "results/susie_coloc_1kg")
LD_FILE_REF <- file.path(FM_DIR, "results/susie_coloc/ld_contamination_clusters.csv")  # read-only ref

cat("============================================================\n")
cat("Combining SuSiE+ABF COLOC results (1KG EUR LD panel)\n")
cat("Input dir :", COLOC_DIR, "\n")
cat("============================================================\n")

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
  cat("  PP.H4.abf > 0.3:", sum(gwas_results$PP.H4.abf > 0.3, na.rm = TRUE), "\n")
  if ("PP.H4.susie" %in% names(gwas_results)) {
    cat("  PP.H4.susie > 0.8:", sum(gwas_results$PP.H4.susie > 0.8, na.rm = TRUE), "\n")
    cat("  PP.H4.susie > 0.5:", sum(gwas_results$PP.H4.susie > 0.5, na.rm = TRUE), "\n")
  }

  fwrite(gwas_results, file.path(gdir, "susie_coloc_combined.csv"))
  all_results[[length(all_results) + 1]] <- gwas_results
}

# Master table
if (length(all_results) > 0) {
  master <- rbindlist(all_results, fill = TRUE)

  # Apply minimum SNP filter (raised from 10 to 100; see threshold audit)
  n_before <- nrow(master)
  master <- master[is.na(n_snps) | n_snps >= MIN_SNPS]
  n_filtered <- n_before - nrow(master)
  cat(sprintf("\n  Filtered %d gene-GWAS pairs with n_snps < %d\n", n_filtered, MIN_SNPS))

  master_file <- file.path(COLOC_DIR, "susie_coloc_all_gwas_1kg.csv")
  fwrite(master, master_file)

  cat("\n============================================================\n")
  cat("Master table:", nrow(master), "gene-GWAS entries\n")
  cat("Written to:", master_file, "\n")

  # ── Sample-overlap-aware GWAS group assignments ──────────────────────────────
  # Identical to v1 (07_combine_susie_coloc.R); only EUR GWAS are present in
  # the 1KG re-run but we keep the full table so unrecognised names still get
  # their own (independent) group via the fallback below.
  # ─────────────────────────────────────────────────────────────────────────────
  gwas_groups <- data.table(
    gwas_name = c(
      "UKBB_ALT", "UKBB_AST", "UKBB_GGT",
      "2019_31311600_NAFLD_EUR",
      "2020_32298765_NAFLD_EUR",
      "2021_34841290_NAFLD_EUR",
      "2023_36280732_NAFLD_deCode_EUR",
      "2023_36280732_NAFLD_Intermountain_EUR",
      "2023_36280732_NAFLD_UKBB_EUR",
      "2021_34128465_PDFF_EUR",
      "2021_34957434_PDFF_EUR",
      "2022_36402844_PDFF_EUR",
      "FinnGen_NAFLD",
      "FinnGen_NASH",
      # cirrhosis/HCC dropped 2026-06-06 (not MASLD-specific):
      # FinnGen_HCC, Ghouse_Cirrhosis, Ghouse_HCC,
      # 2020_32514122_Cirrhosis_EAS, 2020_32514122_HCC_EAS
      "BBJ_ALT", "BBJ_AST", "BBJ_GGT",
      "PanUKBB_AFR_ALT", "PanUKBB_AFR_AST", "PanUKBB_AFR_GGT",
      "PanUKBB_CSA_ALT", "PanUKBB_CSA_AST", "PanUKBB_CSA_GGT"
    ),
    gwas_group = c(
      "UKBB_ALT", "UKBB_AST", "UKBB_GGT",
      "Ghodsian_meta",
      "Anstee_biopsy",
      "Ghodsian_meta",
      "Sveinbjornsson_deCode",
      "Sveinbjornsson_Intermountain",
      "Ghodsian_meta",
      "UKBB_PDFF", "UKBB_PDFF", "UKBB_PDFF",
      "FinnGen_R12",
      "FinnGen_R12",
      # cirrhosis/HCC dropped 2026-06-06 (not MASLD-specific)
      "BBJ_ALT", "BBJ_AST", "BBJ_GGT",
      "PanUKBB_AFR", "PanUKBB_AFR", "PanUKBB_AFR",
      "PanUKBB_CSA", "PanUKBB_CSA", "PanUKBB_CSA"
    )
  )

  master <- merge(master, gwas_groups, by = "gwas_name", all.x = TRUE)
  master[is.na(gwas_group), gwas_group := gwas_name]

  HAS_SUSIE <- "PP.H4.susie" %in% names(master)
  if (HAS_SUSIE) {
    cat("  SuSiE columns detected — will compute parallel SuSiE aggregations\n")
  } else {
    cat("  No SuSiE columns detected — SuSiE aggregation columns will be NA\n")
  }

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
    mutate(
      coloc_best_pp4 = ifelse(is.infinite(coloc_best_pp4), NA, coloc_best_pp4)
    )

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
      mutate(
        coloc_best_susie_pp4 = ifelse(is.infinite(coloc_best_susie_pp4), NA, coloc_best_susie_pp4)
      )

    gene_best <- left_join(gene_best, gene_susie, by = c("gene", "ensembl"))

    n_susie_05 <- sum(gene_best$coloc_n_gwas_susie_h4_05 > 0, na.rm = TRUE)
    n_susie_08 <- sum(gene_best$coloc_n_gwas_susie_h4_08 > 0, na.rm = TRUE)
    n_susie_09 <- sum(gene_best$coloc_n_gwas_susie_h4_09 > 0, na.rm = TRUE)
    n_susie_rep2_group <- sum(gene_best$coloc_n_groups_susie_h4_05 >= 2, na.rm = TRUE)
    cat(sprintf(
      "  SuSiE PP.H4>0.5: %d | >0.8: %d | >0.9: %d | >=2 indep groups: %d\n",
      n_susie_05, n_susie_08, n_susie_09, n_susie_rep2_group
    ))
  } else {
    gene_best$coloc_best_susie_pp4       <- NA_real_
    gene_best$coloc_best_susie_gwas      <- NA_character_
    gene_best$coloc_n_gwas_susie_h4_05   <- NA_integer_
    gene_best$coloc_n_gwas_susie_h4_08   <- NA_integer_
    gene_best$coloc_n_gwas_susie_h4_09   <- NA_integer_
    gene_best$coloc_susie_n_pairs_total  <- NA_integer_
    gene_best$coloc_susie_success_rate   <- NA_real_
    gene_best$coloc_n_groups_susie_h4_05 <- NA_integer_
    gene_best$coloc_n_groups_susie_h4_08 <- NA_integer_
  }

  n_rep2_gwas  <- sum(gene_best$coloc_n_gwas_h4_05  >= 2, na.rm = TRUE)
  n_rep2_group <- sum(gene_best$coloc_n_groups_h4_05 >= 2, na.rm = TRUE)
  cat(sprintf(
    "\n  Independent group replication (PP.H4>0.5): %d genes in >=2 groups (was %d by raw GWAS count)\n",
    n_rep2_group, n_rep2_gwas
  ))
  n_rep2_gwas08  <- sum(gene_best$coloc_n_gwas_h4_08  >= 2, na.rm = TRUE)
  n_rep2_group08 <- sum(gene_best$coloc_n_groups_h4_08 >= 2, na.rm = TRUE)
  cat(sprintf(
    "  Independent group replication (PP.H4>0.8): %d genes in >=2 groups (was %d by raw GWAS count)\n",
    n_rep2_group08, n_rep2_gwas08
  ))

  # ── MHC region flag (chr6:25-35Mb, hg19 coordinates) ──
  MHC_CHR <- 6L; MHC_START <- 25000000L; MHC_END <- 35000000L
  gene_chr <- master[, .(chr = chr[1]), by = .(gene, ensembl)]
  mhc_genes <- master[chr == MHC_CHR & !is.na(top_snp)]
  if (nrow(mhc_genes) > 0) {
    mhc_genes[, snp_pos := as.numeric(sub(".*:", "", top_snp))]
    mhc_ensembl <- unique(mhc_genes[snp_pos >= MHC_START & snp_pos <= MHC_END, ensembl])
  } else {
    mhc_ensembl <- character(0)
  }
  gene_best$is_mhc <- gene_best$ensembl %in% mhc_ensembl
  cat("  MHC-flagged genes:", sum(gene_best$is_mhc), "\n")

  # ── LD contamination cluster annotation (shared annotation; read from v1 dir) ──
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
    n_ld <- sum(!is.na(gene_best$ld_cluster_flag))
    cat("  LD-cluster-flagged genes:", n_ld, "(annotation source:", LD_FILE_REF, ")\n")
  } else {
    gene_best[, ld_cluster_best_gene := NA_character_]
    gene_best[, ld_cluster_flag := NA_character_]
    cat("  LD contamination file not found at", LD_FILE_REF, "— columns added as NA\n")
  }

  gene_file <- file.path(COLOC_DIR, "gene_level_coloc_1kg.csv")
  fwrite(gene_best, gene_file, na = "NA")
  cat("Gene-level table:", nrow(gene_best), "genes\n")
  cat("Written to:", gene_file, "\n")
  cat("============================================================\n")
}
