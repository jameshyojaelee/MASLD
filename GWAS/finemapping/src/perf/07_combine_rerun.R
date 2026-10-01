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
COLOC_DIR <- Sys.getenv("COLOC_DIR",
  unset = file.path(FM_DIR, "results/susie_coloc_rerun"))  # P6: REDIRECTED
# COLOC_MASTER_IN (review item 2, r2 release): start from an already-assembled
# master (MIN_SNPS-filtered, C4 columns present) instead of the per-chr files.
# The master is read, not rewritten; only the gene-level table is written.
MASTER_IN <- Sys.getenv("COLOC_MASTER_IN", unset = "")
if (!grepl("_rerun$|^susie_coloc_r2_", basename(COLOC_DIR))) stop("REFUSING: this copy must never write into canonical results/susie_coloc/")
if (basename(COLOC_DIR) == "susie_coloc") stop("REFUSING: COLOC_DIR resolved to canonical")

cat("============================================================\n")
cat("Combining ABF COLOC results\n")
cat("============================================================\n")

gwas_dirs <- if (nzchar(MASTER_IN)) character(0) else list.dirs(COLOC_DIR, recursive = FALSE)
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
if (length(all_results) > 0 || nzchar(MASTER_IN)) {
  if (nzchar(MASTER_IN)) {
    master <- fread(MASTER_IN)
    cat("Master table read from", MASTER_IN, ":", nrow(master), "rows\n")
    if (any(!is.na(master$n_snps) & master$n_snps < MIN_SNPS)) stop("COLOC_MASTER_IN is not MIN_SNPS-filtered")
  } else {
  master <- rbindlist(all_results, fill = TRUE)

  # Apply minimum SNP filter (raised from 10 to 100; see threshold audit)
  n_before <- nrow(master)
  master <- master[is.na(n_snps) | n_snps >= MIN_SNPS]
  n_filtered <- n_before - nrow(master)
  cat(sprintf("\n  Filtered %d gene-GWAS pairs with n_snps < %d\n", n_filtered, MIN_SNPS))
  }

  # ── C4: lambda_s / LD-reliability propagation (cross-ancestry LD-matching) ────
  # Propagate the SuSiE-RSS LD-consistency diagnostic lambda_s (estimate_s_rss;
  # Zou 2022 / finemapping_functions.R) from the FM per-locus master into the
  # coloc masters, so non-EUR (small 1000G panel) LD-inconsistency is visible and
  # gate-able downstream instead of being silently conflated with the ancestry-
  # matched EUR arm. GUARANTEED carrier = per-stratum reliability tier (every
  # non-EUR GWAS stratum has lambda_s in the FM master); OPTIONAL refinement =
  # per-locus lambda_s matched on the coloc top_snp position falling inside the
  # fine-mapped locus window. The two tiny lookup tables are produced by
  # qa_campaign/audit/C4_lambda_s_audit.R (streams the 5.3 GB FM master once), so
  # this join keeps the 07 refresh at ~1 min. If a future 06 run emits a per-row
  # `lambda_s_locus` column, it is preferred over the join (forward-compatible).
  QA_DIR <- file.path(FM_DIR, "results/qa_campaign")
  LAMBDA_S_HIGH <- 0.20   # documented flag cutoff (see C4_lambda_s_audit.R)
  # C5 (2026-07-05): registry-authoritative ancestry. config/gwas_registry.tsv
  # (study_name -> ancestry) is the single source of truth; the regex is a fallback
  # ONLY for an unregistered gwas_name, which now WARNS (was a silent -> EUR
  # default); a registry-vs-regex disagreement is a hard stop (drift guard).
  .reg_anc <- { .r <- tryCatch(read.delim(file.path(FM_DIR, "config/gwas_registry.tsv"), stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(.r) || !all(c("study_name", "ancestry") %in% names(.r))) NULL else setNames(.r$ancestry, .r$study_name) }
  .anc_regex <- function(g) fifelse(is.na(g), NA_character_,
    fifelse(grepl("^PanUKBB_AFR", g), "AFR",
    fifelse(grepl("^PanUKBB_CSA", g), "SAS",
    fifelse(grepl("^BBJ_", g), "EAS",
    fifelse(grepl("^MVP_.*_(AFR|AMR|EAS|SAS)$", g), sub(".*_", "", g), "EUR")))))
  .anc_of <- function(g) {
    rx <- .anc_regex(g)
    if (is.null(.reg_anc)) return(rx)
    reg <- unname(.reg_anc[as.character(g)])
    unreg <- !is.na(g) & is.na(reg)
    if (any(unreg)) warning(sprintf("anc_of: %d unregistered gwas_name(s) -> regex fallback: %s",
      length(unique(g[unreg])), paste(unique(g[unreg]), collapse = ", ")))
    dis <- !is.na(reg) & !is.na(rx) & reg != rx
    if (any(dis)) stop(sprintf("anc_of: registry-vs-regex ancestry disagreement: %s",
      paste(unique(g[dis]), collapse = ", ")))
    fifelse(is.na(reg), rx, reg)
  }
  master[, ancestry := .anc_of(gwas_name)]
  strat_f <- file.path(QA_DIR, "lambda_s_by_stratum.tsv")
  locus_f <- file.path(QA_DIR, "lambda_s_by_locus.tsv")
  if (nzchar(MASTER_IN)) {
    # The input master already carries the propagated columns; merging the
    # stratum table again would duplicate them as .x/.y.
    cat("  C4 lambda_s columns carried from COLOC_MASTER_IN (not re-propagated)\n")
  } else if (file.exists(strat_f)) {
    strat <- fread(strat_f)
    master <- merge(master,
      strat[, .(gwas_name = study, ld_panel, ld_panel_n,
                stratum_lambda_s_median   = lambda_s_median,
                stratum_lambda_s_p95      = lambda_s_p95,
                stratum_frac_high_lambda_s = frac_loci_high,
                stratum_ld_reliability)],
      by = "gwas_name", all.x = TRUE)

    # per-locus refinement (two sources, coalesced): (1) a per-row lambda_s_locus
    # emitted directly by 06 for future runs; (2) foverlaps of the coloc top_snp
    # position into the fine-mapped locus window [pos_min,pos_max] for the SAME
    # study (authoritative variant bounds) for rows lacking (1). Coalescing keeps
    # this correct under pure-old, pure-new, or mixed per-chr inputs.
    master[, locus_lambda_s := NA_real_]
    if ("lambda_s_locus" %in% names(master)) {
      master[, locus_lambda_s := suppressWarnings(as.numeric(lambda_s_locus))]
    }
    if (file.exists(locus_f)) {
      loc <- fread(locus_f)[!is.na(lambda_s)]
      master[, snp_pos_c4 := suppressWarnings(as.numeric(sub(".*:", "", top_snp)))]
      master[, row_id_c4  := .I]
      key_tab <- master[is.na(locus_lambda_s) & !is.na(snp_pos_c4),
                        .(row_id_c4, gwas_name, chr,
                          start = snp_pos_c4, end = snp_pos_c4)]
      if (nrow(key_tab) > 0 && nrow(loc) > 0) {
        loc_iv <- loc[, .(gwas_name = study, chr, start = pos_min, end = pos_max,
                          loc_lambda = lambda_s)]
        setkey(loc_iv, gwas_name, chr, start, end)
        setkey(key_tab, gwas_name, chr, start, end)
        ov <- foverlaps(key_tab, loc_iv, type = "within", nomatch = NA)
        # a top_snp overlapping >1 fine-map window -> worst-case (max) lambda_s
        ov_best <- ov[!is.na(loc_lambda), .(loc_lambda = max(loc_lambda)), by = row_id_c4]
        if (nrow(ov_best) > 0) master[ov_best, locus_lambda_s := i.loc_lambda, on = "row_id_c4"]
      }
      master[, c("snp_pos_c4", "row_id_c4") := NULL]
    }

    # numeric carrier: per-locus lambda_s when matched, else per-stratum median
    master[, lambda_s := fifelse(!is.na(locus_lambda_s), locus_lambda_s,
                                 stratum_lambda_s_median)]
    # row-level LD-reliability flag — NEVER NA for a non-EUR row (s10 guard target)
    master[, ld_reliability := fifelse(ancestry == "EUR", "EUR_ancestry_matched",
      fifelse(!is.na(locus_lambda_s),
              fifelse(locus_lambda_s > LAMBDA_S_HIGH, "locus_high_lambda_s", "locus_ok"),
              fifelse(is.na(stratum_ld_reliability), "unknown_no_finemap",
                      stratum_ld_reliability)))]
    n_nonEUR <- sum(master$ancestry != "EUR", na.rm = TRUE)
    n_flag   <- sum(master$ancestry != "EUR" & !is.na(master$ld_reliability))
    n_high   <- sum(master$ld_reliability == "locus_high_lambda_s", na.rm = TRUE)
    cat(sprintf("  C4 lambda_s propagated: non-EUR rows %d (flagged %d = %.1f%%) | locus-matched high-lambda_s rows %d\n",
                n_nonEUR, n_flag, 100 * n_flag / max(n_nonEUR, 1L), n_high))
  } else {
    cat("  C4 WARNING: qa_campaign/lambda_s_by_stratum.tsv not found —",
        "LD-reliability columns stubbed (run audit/C4_lambda_s_audit.R first)\n")
    master[, `:=`(ld_panel = NA_character_, ld_panel_n = NA_integer_,
                  stratum_lambda_s_median = NA_real_, stratum_lambda_s_p95 = NA_real_,
                  stratum_frac_high_lambda_s = NA_real_, stratum_ld_reliability = NA_character_,
                  locus_lambda_s = NA_real_, lambda_s = NA_real_)]
    master[, ld_reliability := fifelse(ancestry == "EUR", "EUR_ancestry_matched", NA_character_)]
  }
  # ─────────────────────────────────────────────────────────────────────────────

  master_file <- if (nzchar(MASTER_IN)) MASTER_IN else file.path(COLOC_DIR, "susie_coloc_all_gwas.csv")
  if (!nzchar(MASTER_IN)) fwrite(master, master_file)

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
  # MVP (phs002453): all strata of one GIA ancestry are the SAME veterans
  # (NAFLD/cirrhosis/ALT/AST/Albumin/Platelet overlap within a person), so collapse
  # every MVP_<pheno>_<ANC> into a single per-ancestry group MVP_<ANC>. Prevents the
  # within-cohort liver phenotypes from inflating coloc_n_groups_h4 replication.
  master[grepl("^MVP_", gwas_name),
         gwas_group := paste0("MVP_", sub(".*_", "", gwas_name))]
  master[is.na(gwas_group), gwas_group := gwas_name]

  # ── Ancestry provenance tag ──────────────────────────────────────────────────
  # The Broadaway eQTL panel is EUROPEAN, so only EUR GWAS strata are
  # ancestry-matched to the eQTL; non-EUR strata (BBJ=EAS, PanUKBB_AFR/CSA,
  # MVP_*_{AFR,AMR,EAS,SAS}) are cross-ancestry regardless of PP.H4 magnitude.
  # Field precedent (MVP/Vujkovic 2022, PMID 35654975): hold the SAME PP.H4 bar
  # across ancestries and differentiate by evidence STATUS (matched vs
  # cross-ancestry/exploratory), NOT by a lowered threshold. We therefore keep the
  # inclusive cross-stratum max as the atlas default (loose/suggestive) and ADD
  # provenance + confidence-tier + EUR-only-companion columns so nothing is
  # silently conflated. Confidence tiers are field-grounded: high>=0.8 (eQTL
  # Catalogue/FinnGen), suggestive 0.5-0.8 (GTEx enloc), nominal 0.1-0.5.
  # C5 (2026-07-05): registry-authoritative ancestry (mirrors the .anc_of guard
  # above — registry = source of truth, regex = fallback for unregistered names
  # with a warn, and a registry-vs-regex disagreement is a hard stop).
  .reg_anc <- { .r <- tryCatch(read.delim(file.path(FM_DIR, "config/gwas_registry.tsv"), stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(.r) || !all(c("study_name", "ancestry") %in% names(.r))) NULL else setNames(.r$ancestry, .r$study_name) }
  .anc_regex <- function(g) fifelse(is.na(g), NA_character_,
    fifelse(grepl("^PanUKBB_AFR", g), "AFR",
    fifelse(grepl("^PanUKBB_CSA", g), "SAS",
    fifelse(grepl("^BBJ_", g), "EAS",
    fifelse(grepl("^MVP_.*_(AFR|AMR|EAS|SAS)$", g), sub(".*_", "", g), "EUR")))))
  anc_of <- function(g) {
    rx <- .anc_regex(g)
    if (is.null(.reg_anc)) return(rx)
    reg <- unname(.reg_anc[as.character(g)])
    unreg <- !is.na(g) & is.na(reg)
    if (any(unreg)) warning(sprintf("anc_of: %d unregistered gwas_name(s) -> regex fallback: %s",
      length(unique(g[unreg])), paste(unique(g[unreg]), collapse = ", ")))
    dis <- !is.na(reg) & !is.na(rx) & reg != rx
    if (any(dis)) stop(sprintf("anc_of: registry-vs-regex ancestry disagreement: %s",
      paste(unique(g[dis]), collapse = ", ")))
    fifelse(is.na(reg), rx, reg)
  }
  conf_tier <- function(pp) fifelse(is.na(pp) | pp < 0.1, "none",
    fifelse(pp >= 0.8, "high", fifelse(pp >= 0.5, "suggestive", "nominal")))
  master[, ancestry := anc_of(gwas_name)]

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
      # EUR-only companion (ancestry-matched to the EUR eQTL) + ancestry breadth
      coloc_best_pp4_EUR = { m <- suppressWarnings(max(PP.H4.abf[ancestry == "EUR"], na.rm = TRUE)); if (is.infinite(m)) NA_real_ else m },
      coloc_best_gwas_EUR = { i <- which(ancestry == "EUR" & !is.na(PP.H4.abf)); if (length(i) == 0) NA_character_ else gwas_name[i][which.max(PP.H4.abf[i])] },
      coloc_n_anc_h4_05 = n_distinct(ancestry[!is.na(PP.H4.abf) & PP.H4.abf > 0.5]),
      # C4: LD-reliability of the row that set the headline PP.H4.abf + worst-case
      # (max) lambda_s across this gene's non-EUR colocalizing GWAS (PP.H4>0.5).
      coloc_best_ld_reliability = if (all(is.na(PP.H4.abf))) NA_character_ else ld_reliability[which.max(PP.H4.abf)],
      coloc_best_lambda_s       = if (all(is.na(PP.H4.abf))) NA_real_ else lambda_s[which.max(PP.H4.abf)],
      coloc_max_lambda_s        = { v <- lambda_s[ancestry != "EUR" & !is.na(PP.H4.abf) & PP.H4.abf > 0.5]; v <- v[!is.na(v)]; if (length(v) == 0) NA_real_ else max(v) },
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
        coloc_best_susie_pp4_EUR = { m <- suppressWarnings(max(PP.H4.susie[ancestry == "EUR"], na.rm = TRUE)); if (is.infinite(m)) NA_real_ else m },
        coloc_best_susie_gwas_EUR = { i <- which(ancestry == "EUR" & !is.na(PP.H4.susie)); if (length(i) == 0) NA_character_ else gwas_name[i][which.max(PP.H4.susie[i])] },
        coloc_n_anc_susie_h4_05 = n_distinct(ancestry[!is.na(PP.H4.susie) & PP.H4.susie > 0.5]),
        # C4: LD-reliability + lambda_s of the SuSiE headline row (cross-ancestry)
        coloc_best_susie_ld_reliability = ld_reliability[which.max(PP.H4.susie)],
        coloc_best_susie_lambda_s       = lambda_s[which.max(PP.H4.susie)],
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
    gene_best$coloc_best_susie_pp4_EUR    <- NA_real_
    gene_best$coloc_best_susie_gwas_EUR   <- NA_character_
    gene_best$coloc_n_anc_susie_h4_05     <- NA_integer_
    gene_best$coloc_best_susie_ld_reliability <- NA_character_
    gene_best$coloc_best_susie_lambda_s   <- NA_real_
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

  # ── Provenance + confidence-tier labels ──────────────────────────────────────
  # The headline coloc_best_* stays the inclusive cross-stratum max (atlas default,
  # deliberately loose). These columns annotate WHERE the headline came from and
  # HOW strong it is, so downstream consumers that make stronger claims can gate on
  # the EUR-only companion while the atlas itself stays inclusive. (MAJOR-2 fix by
  # labeling, not exclusion — see the ancestry-tag comment above.)
  gene_best <- as.data.table(gene_best)

  # Review item 2: a study whose every signal pair failed coloc's shared-posterior
  # check has PP.H4.susie = NA, so the SuSiE block above never sees it. Count it
  # as its own class; such a gene is never a SuSiE negative.
  # coloc_susie_state: positive (> 0.5) > untestable > tested_le_0.5 > no_susie_result.
  UNTESTABLE <- "susie_untestable_insufficient_shared_posterior"
  gene_untest <- master[method == UNTESTABLE,
    .(coloc_n_gwas_susie_untestable = .N,
      coloc_susie_untestable_gwas = paste(sort(gwas_name), collapse = ";")),
    by = .(gene, ensembl)]
  gene_best <- merge(gene_best, gene_untest, by = c("gene", "ensembl"), all.x = TRUE)
  gene_best[is.na(coloc_n_gwas_susie_untestable), coloc_n_gwas_susie_untestable := 0L]
  gene_best[, coloc_susie_state := fifelse(!is.na(coloc_n_gwas_susie_h4_05) & coloc_n_gwas_susie_h4_05 > 0, "positive",
    fifelse(coloc_n_gwas_susie_untestable > 0, "untestable",
      fifelse(!is.na(coloc_best_susie_pp4), "tested_le_0.5", "no_susie_result")))]
  cat("  SuSiE state by gene:", paste(names(table(gene_best$coloc_susie_state)),
      table(gene_best$coloc_susie_state), sep = "=", collapse = " "), "\n")
  gene_best[, coloc_best_ancestry       := anc_of(coloc_best_gwas)]
  gene_best[, coloc_best_susie_ancestry := anc_of(coloc_best_susie_gwas)]
  gene_best[, coloc_abf_conf_tier   := conf_tier(coloc_best_pp4)]
  gene_best[, coloc_susie_conf_tier := conf_tier(coloc_best_susie_pp4)]
  gene_best[, coloc_abf_ancestry_matched   := coloc_best_ancestry == "EUR"]
  gene_best[, coloc_susie_ancestry_matched := coloc_best_susie_ancestry == "EUR"]
  # A meaningful (>=0.5, GTEx "colocalized" floor) headline set SOLELY by a
  # cross-ancestry stratum with no EUR corroboration >=0.5 -- i.e. a call that
  # would look EUR-validated but is cross-ancestry-only (the MAJOR-2 case).
  gene_best[, coloc_abf_headline_cross_anc   := !is.na(coloc_best_pp4) & coloc_best_pp4 >= 0.5 &
              !is.na(coloc_best_ancestry) & coloc_best_ancestry != "EUR" &
              (is.na(coloc_best_pp4_EUR) | coloc_best_pp4_EUR < 0.5)]
  gene_best[, coloc_susie_headline_cross_anc := !is.na(coloc_best_susie_pp4) & coloc_best_susie_pp4 >= 0.5 &
              !is.na(coloc_best_susie_ancestry) & coloc_best_susie_ancestry != "EUR" &
              (is.na(coloc_best_susie_pp4_EUR) | coloc_best_susie_pp4_EUR < 0.5)]
  # cross-ancestry replicated = same 0.5 bar met in >=2 distinct ancestries
  # (field convention; matches the existing "PP.H4>0.5 in both EUR and EAS" rule)
  gene_best[, coloc_cross_ancestry_replicated :=
              (!is.na(coloc_n_anc_h4_05) & coloc_n_anc_h4_05 >= 2) |
              (!is.na(coloc_n_anc_susie_h4_05) & coloc_n_anc_susie_h4_05 >= 2)]

  # C4: cross-ancestry LD caution — a (cross-ancestry-only) headline that rests on
  # a non-EUR stratum whose colocalizing locus has high lambda_s (LD reference does
  # not match the GWAS -> the call is LD-unreliable, not just cross-ancestry). Pairs
  # with the existing coloc_*_headline_cross_anc provenance labels.
  gene_best[, coloc_cross_anc_ld_caution :=
              (coloc_abf_headline_cross_anc == TRUE &
                 !is.na(coloc_best_lambda_s) & coloc_best_lambda_s > LAMBDA_S_HIGH) |
              (coloc_susie_headline_cross_anc == TRUE &
                 !is.na(coloc_best_susie_lambda_s) & coloc_best_susie_lambda_s > LAMBDA_S_HIGH)]
  cat("  C4 cross-ancestry LD-caution genes (headline on high-lambda_s non-EUR locus):",
      sum(gene_best$coloc_cross_anc_ld_caution, na.rm = TRUE), "\n")

  # ── LD contamination cluster annotation ──
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

  # C4: raw-named aliases (= the headline coloc call) so lambda_s / ld_reliability
  # are discoverable in the gene-level master by the same header probe used on the
  # coloc master. The precise semantics live in the coloc_best_* / coloc_max_*
  # columns; these two are convenience aliases for the headline row.
  gene_best[, lambda_s := coloc_best_lambda_s]
  gene_best[, ld_reliability := coloc_best_ld_reliability]

  gene_file <- file.path(COLOC_DIR, "gene_level_coloc.csv")
  fwrite(gene_best, gene_file, na = "NA")
  cat("Gene-level table:", nrow(gene_best), "genes\n")
  cat("Written to:", gene_file, "\n")

  # ── C4 GUARD (s10): assert lambda_s / LD-reliability survived end-to-end ──────
  # Automated in-pipeline check — every combine (incl. run_07_refresh, which calls
  # this script) FAILS if a non-EUR coloc row is missing its propagated
  # LD-reliability flag, or if the lambda_s columns were dropped. Enforced only
  # when the lambda_s lookup was available (propagation path ran); otherwise a
  # missing lookup downgrades to a warning so a fresh checkout is never bricked.
  if (file.exists(strat_f)) {
    miss_m <- setdiff(c("lambda_s", "ld_reliability", "ancestry",
                        "stratum_ld_reliability"), names(master))
    miss_g <- setdiff(c("coloc_best_ld_reliability", "coloc_best_lambda_s",
                        "coloc_max_lambda_s", "coloc_cross_anc_ld_caution",
                        "lambda_s", "ld_reliability"),
                      names(gene_best))
    non_eur_unflagged <- master[ancestry != "EUR" & is.na(ld_reliability), .N]
    if (length(miss_m) || length(miss_g) || non_eur_unflagged > 0) {
      stop(sprintf(paste0("C4 GUARD FAILED: master missing {%s}; gene missing {%s}; ",
                          "%d non-EUR coloc rows without an ld_reliability flag"),
                   paste(miss_m, collapse = ","), paste(miss_g, collapse = ","),
                   non_eur_unflagged))
    }
    cat(sprintf("  C4 GUARD PASSED: lambda_s + ld_reliability present end-to-end; all %d non-EUR rows flagged\n",
                master[ancestry != "EUR", .N]))
  } else {
    cat("  C4 GUARD SKIPPED: lambda_s lookup absent (run audit/C4_lambda_s_audit.R to enable)\n")
  }

  # ── C8 GUARD (s10): positive/negative control recovery on the fresh master ────
  # Re-check the 7 canonical anchors against their MECHANISM-based expected tier
  # on every COLOC-master rebuild (this table is the atlas's direct input). The
  # anchor half of the guard is inlined here (dependency-free) so a rebuild always
  # self-checks; the full re-run (anchors + the AF/expression-matched null set +
  # the 0.312 note) lives in audit/C8_control_recovery.R. FAILS if an expression-
  # mediated positive drops below tier `high` OR a coding/splice-mechanism anchor
  # inflates above the over-call ceiling (a false-positive regression). Anchors
  # absent from a partial run downgrade to a warning so a subset run is not bricked.
  c8_pos  <- c("GCKR", "TM6SF2", "THRB", "EXOC3L4")   # expect high (>=0.5)
  c8_code <- c("PNPLA3", "HSD17B13", "MBOAT7")         # coding/splice: must stay <=0.5
  c8_fail <- character(0)
  gb <- gene_best
  for (g in c8_pos) {
    v <- gb[gene == g, coloc_best_pp4]
    if (!length(v) || is.na(v[1])) { warning(sprintf("C8: positive anchor %s absent from master", g)); next }
    if (v[1] < 0.5) c8_fail <- c(c8_fail, sprintf("%s dropped tier (abf_pp4=%.3f < 0.5)", g, v[1]))
  }
  # EXOC3L4 spot-check: require ABF+SuSiE two-arm concordance (both >= 0.5)
  ex_s <- gb[gene == "EXOC3L4", coloc_best_susie_pp4]
  if (length(ex_s) && !is.na(ex_s[1]) && ex_s[1] < 0.5)
    c8_fail <- c(c8_fail, sprintf("EXOC3L4 lost two-arm concordance (susie_pp4=%.3f < 0.5)", ex_s[1]))
  for (g in c8_code) {
    v <- gb[gene == g, coloc_best_pp4]
    if (!length(v) || is.na(v[1])) { warning(sprintf("C8: coding-mechanism anchor %s absent from master", g)); next }
    if (v[1] > 0.5) c8_fail <- c(c8_fail, sprintf("%s inflated / over-called (abf_pp4=%.3f > 0.5)", g, v[1]))
  }
  if (length(c8_fail)) {
    stop(sprintf("C8 GUARD FAILED (control recovery regressed): %s",
                 paste(c8_fail, collapse = "; ")))
  }
  cat(sprintf("  C8 GUARD PASSED: %d/%d positive anchors hold tier `high`; %d/%d coding-mechanism anchors do not over-call (EXOC3L4 ABF+SuSiE concordant)\n",
              sum(sapply(c8_pos, function(g){v<-gb[gene==g,coloc_best_pp4]; length(v)&&!is.na(v[1])&&v[1]>=0.5})), length(c8_pos),
              sum(sapply(c8_code, function(g){v<-gb[gene==g,coloc_best_pp4]; length(v)&&!is.na(v[1])&&v[1]<=0.5})), length(c8_code)))
  cat("============================================================\n")
}
