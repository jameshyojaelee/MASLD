#!/usr/bin/env Rscript
# 07_combine_susie_coloc.R
# Combines per-chromosome ABF COLOC results across all GWAS
# Usage: Rscript 07_combine_susie_coloc.R

library(data.table)
library(dplyr)

MIN_SNPS <- 100L  # Post-hoc filter: exclude gene-GWAS pairs with < 100 overlapping SNPs

BASE_DIR <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FM_DIR <- file.path(BASE_DIR, "GWAS/finemapping")
COLOC_DIR <- file.path(FM_DIR, "results/susie_coloc")

cat("============================================================\n")
cat("Combining ABF COLOC results\n")
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
  cat("  PP.H4 > 0.8:", sum(gwas_results$PP.H4.abf > 0.8, na.rm = TRUE), "\n")
  cat("  PP.H4 > 0.5:", sum(gwas_results$PP.H4.abf > 0.5, na.rm = TRUE), "\n")
  cat("  PP.H4 > 0.3:", sum(gwas_results$PP.H4.abf > 0.3, na.rm = TRUE), "\n")

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

  master_file <- file.path(COLOC_DIR, "susie_coloc_all_gwas.csv")
  fwrite(master, master_file)

  cat("\n============================================================\n")
  cat("Master table:", nrow(master), "gene-GWAS entries\n")
  cat("Written to:", master_file, "\n")

  # ── Sample-overlap-aware GWAS group assignments ──────────────────────────────
  # Several studies share participants (e.g., three UKBB liver-enzyme GWAS share
  # ~344K individuals; three PDFF GWAS share ~33K UKBB MRI subjects; Ghodsian
  # subsumes the eMERGE cohort and shares UKBB controls with Sveinbjornsson_UKBB).
  # Counting raw GWAS inflates the apparent replication signal.
  #
  # Strategy: assign each GWAS to a gwas_group reflecting independent participant
  # pools.  coloc_n_groups_h4_05/08 count UNIQUE groups with PP.H4 above the
  # threshold and should be preferred over coloc_n_gwas_h4_05/08 for any claim
  # of "replicated across independent studies".
  #
  # coloc_n_gwas_h4_05/08 are retained unchanged for backward compatibility.
  # ─────────────────────────────────────────────────────────────────────────────
  gwas_groups <- data.table(
    gwas_name = c(
      # UKBB liver enzymes: same ~344K participants but genuinely different traits
      "UKBB_ALT", "UKBB_AST", "UKBB_GGT",
      # NAFLD case-control
      "2019_31311600_NAFLD_EUR",               # eMERGE — subsumed by Ghodsian
      "2020_32298765_NAFLD_EUR",               # Anstee biopsy — independent
      "2021_34841290_NAFLD_EUR",               # Ghodsian meta (eMERGE+UKBB+Estonian+FinnGen DF4)
      "2023_36280732_NAFLD_deCode_EUR",        # Sveinbjornsson deCODE — independent
      "2023_36280732_NAFLD_Intermountain_EUR", # Sveinbjornsson Intermountain — independent
      "2023_36280732_NAFLD_UKBB_EUR",          # Sveinbjornsson UKBB — overlaps Ghodsian UKBB arm
      # PDFF: all ~33K UKBB MRI participants, different phenotyping pipelines
      "2021_34128465_PDFF_EUR",  # Liu abdominal MRI
      "2021_34957434_PDFF_EUR",  # Haas ML-derived
      "2022_36402844_PDFF_EUR",  # van der Meer whole-body MRI
      # FinnGen R12: restored 2026-04-09; overlaps Ghodsian DF4 arm (651/4614 cases = 14%)
      "FinnGen_NAFLD",            # FinnGen R12 NAFLD (4,614 cases)
      "FinnGen_NASH",             # FinnGen R12 NASH (1,823 cases)
      # Cirrhosis/HCC GWAS (FinnGen_HCC, Ghouse_Cirrhosis/HCC, Ishigaki EAS) dropped
      # 2026-06-06: etiology-mixed endpoints excluded by Broadaway; not MASLD-specific.
      # BBJ EAS liver enzymes: same ~160K BBJ participants, different traits
      "BBJ_ALT", "BBJ_AST", "BBJ_GGT",
      # Pan-UKBB AFR liver enzymes: same ~6.6K AFR participants
      "PanUKBB_AFR_ALT", "PanUKBB_AFR_AST", "PanUKBB_AFR_GGT",
      # Pan-UKBB CSA liver enzymes: same ~8.9K CSA participants
      "PanUKBB_CSA_ALT", "PanUKBB_CSA_AST", "PanUKBB_CSA_GGT"
    ),
    gwas_group = c(
      "UKBB_ALT", "UKBB_AST", "UKBB_GGT",
      "Ghodsian_meta",           # eMERGE subsumed by Ghodsian
      "Anstee_biopsy",           # independent biopsy cohort
      "Ghodsian_meta",           # Ghodsian is the primary large NAFLD GWAS
      "Sveinbjornsson_deCode",   # independent Icelandic cohort
      "Sveinbjornsson_Intermountain", # independent US cohort
      "Ghodsian_meta",           # UKBB arm overlaps Ghodsian UKBB participants
      "UKBB_PDFF", "UKBB_PDFF", "UKBB_PDFF",
      "FinnGen_R12",             # FinnGen R12 overlaps Ghodsian DF4 — own group
      "FinnGen_R12",             # Same FinnGen R12 cohort, different phenotype
      "BBJ_ALT", "BBJ_AST", "BBJ_GGT",
      "PanUKBB_AFR", "PanUKBB_AFR", "PanUKBB_AFR",
      "PanUKBB_CSA", "PanUKBB_CSA", "PanUKBB_CSA"
    )
  )

  # Merge gwas_group into master; unrecognised GWAS names get their own group
  # (conservative: treats unknown studies as independent)
  master <- merge(master, gwas_groups, by = "gwas_name", all.x = TRUE)
  master[is.na(gwas_group), gwas_group := gwas_name]

  # Guard: detect whether SuSiE columns are present in this run
  HAS_SUSIE <- "PP.H4.susie" %in% names(master)
  if (HAS_SUSIE) {
    cat("  SuSiE columns detected — will compute parallel SuSiE aggregations\n")
  } else {
    cat("  No SuSiE columns detected — SuSiE aggregation columns will be NA\n")
  }

  # Best PP.H4 per gene across all GWAS (for atlas integration)
  gene_best <- master %>%
    group_by(gene, ensembl) %>%
    summarise(
      coloc_best_pp4 = max(PP.H4.abf, na.rm = TRUE),
      coloc_best_gwas = if (all(is.na(PP.H4.abf))) NA_character_ else gwas_name[which.max(PP.H4.abf)],
      coloc_n_gwas_h4_05 = sum(PP.H4.abf > 0.5, na.rm = TRUE),
      coloc_n_gwas_h4_08 = sum(PP.H4.abf > 0.8, na.rm = TRUE),
      coloc_n_gwas_tested = n(),
      # Group-level counts (sample-overlap-corrected)
      coloc_n_groups_h4_05 = n_distinct(gwas_group[!is.na(PP.H4.abf) & PP.H4.abf > 0.5]),
      coloc_n_groups_h4_08 = n_distinct(gwas_group[!is.na(PP.H4.abf) & PP.H4.abf > 0.8]),
      .groups = "drop"
    ) %>%
    mutate(
      coloc_best_pp4 = ifelse(is.infinite(coloc_best_pp4), NA, coloc_best_pp4)
    )

  # ── SuSiE-specific per-gene aggregation ──────────────────────────────────────
  # Mirrors the ABF aggregation above, including sample-overlap correction via
  # gwas_group.  Only executed when SuSiE columns are present (HAS_SUSIE = TRUE).
  # When absent, stub columns of NA are appended so gene_level_coloc.csv always
  # has the same schema.
  # ─────────────────────────────────────────────────────────────────────────────
  if (HAS_SUSIE) {
    master_susie <- as.data.frame(master)  # dplyr-compatible view
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
        # Group-level counts (same sample-overlap correction as ABF block above)
        coloc_n_groups_susie_h4_05 = n_distinct(gwas_group[!is.na(PP.H4.susie) & PP.H4.susie > 0.5]),
        coloc_n_groups_susie_h4_08 = n_distinct(gwas_group[!is.na(PP.H4.susie) & PP.H4.susie > 0.8]),
        .groups = "drop"
      ) %>%
      mutate(
        coloc_best_susie_pp4 = ifelse(is.infinite(coloc_best_susie_pp4), NA, coloc_best_susie_pp4)
      )

    gene_best <- left_join(gene_best, gene_susie, by = c("gene", "ensembl"))

    # Summary counts
    n_susie_05 <- sum(gene_best$coloc_n_gwas_susie_h4_05 > 0, na.rm = TRUE)
    n_susie_08 <- sum(gene_best$coloc_n_gwas_susie_h4_08 > 0, na.rm = TRUE)
    n_susie_rep2_group <- sum(gene_best$coloc_n_groups_susie_h4_05 >= 2, na.rm = TRUE)
    cat(sprintf(
      "  SuSiE PP.H4>0.5: %d genes | PP.H4>0.8: %d genes | >=2 independent groups: %d genes\n",
      n_susie_05, n_susie_08, n_susie_rep2_group
    ))
  } else {
    # Stub columns so downstream scripts always see a consistent schema
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
  # ─────────────────────────────────────────────────────────────────────────────

  # Summary: compare raw-GWAS vs group-based replication counts
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
  # Get per-gene chr from master table (use first occurrence)
  gene_chr <- master[, .(chr = chr[1]), by = .(gene, ensembl)]
  # Mark MHC: need position info — use top_snp to extract position
  mhc_genes <- master[chr == MHC_CHR & !is.na(top_snp)]
  if (nrow(mhc_genes) > 0) {
    mhc_genes[, snp_pos := as.numeric(sub(".*:", "", top_snp))]
    mhc_ensembl <- unique(mhc_genes[snp_pos >= MHC_START & snp_pos <= MHC_END, ensembl])
  } else {
    mhc_ensembl <- character(0)
  }
  gene_best$is_mhc <- gene_best$ensembl %in% mhc_ensembl
  cat("  MHC-flagged genes:", sum(gene_best$is_mhc), "\n")

  # ── LD contamination cluster annotation ──
  gene_best <- as.data.table(gene_best)  # ensure data.table for consistent merge
  ld_file <- file.path(COLOC_DIR, "ld_contamination_clusters.csv")
  if (file.exists(ld_file)) {
    ld_clusters <- fread(ld_file)
    # Expand cluster genes to a flat lookup table
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
    ld_dt <- ld_dt[!duplicated(gene)]  # keep first occurrence per gene
    gene_best <- merge(gene_best, ld_dt, by = "gene", all.x = TRUE)
    n_ld <- sum(!is.na(gene_best$ld_cluster_flag))
    cat("  LD-cluster-flagged genes:", n_ld, "\n")
  } else {
    gene_best[, ld_cluster_best_gene := NA_character_]
    gene_best[, ld_cluster_flag := NA_character_]
    cat("  LD contamination file not found — columns added as NA\n")
  }

  gene_file <- file.path(COLOC_DIR, "gene_level_coloc.csv")
  fwrite(gene_best, gene_file, na = "NA")
  cat("Gene-level table:", nrow(gene_best), "genes\n")
  cat("Written to:", gene_file, "\n")
  cat("============================================================\n")
}
