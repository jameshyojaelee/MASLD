#!/usr/bin/env Rscript
# 08_gwas_provenance_table.R
# ------------------------------------------------------------------------------
# Build the per-GWAS provenance + analysis-coverage supplementary table for the
# paper: every GWAS dataset we hold, its metadata (ancestry / trait / N / case-
# control / LD panel / genome build), and the analyses performed on it
# (fine-mapping: SuSiE + CARMA + SuSiEx/MESuSiE where multi-ancestry; colocalization:
# coloc.abf and SuSiE-coloc, with PP.H4 hit counts at 0.5/0.8/0.9).
#
# Reads ONLY on-disk artifacts so it is re-runnable and self-refreshing:
#   - config/gwas_registry.tsv                (metadata, 1 row per GWAS)
#   - results/combined_finemapping.csv        (fine-mapping, study column)
#   - results/susie_coloc/<study>/susie_coloc_chr*.csv  (per-GWAS COLOC)
#
# Outputs:
#   - results/gwas_analysis_provenance.csv    (machine-readable, all GWAS)
#   - results/gwas_analysis_provenance.md     (paper-ready markdown table)
#
# Usage: Rscript src/08_gwas_provenance_table.R
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({ library(data.table) })

BASE_DIR <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FM_DIR <- file.path(BASE_DIR, "GWAS/finemapping")
setwd(FM_DIR)

PP4_HI <- 0.5   # primary "colocalized" threshold (matches manuscript canonical)

# ── 1. Registry (one row per GWAS) ───────────────────────────────────────────
reg <- fread("config/gwas_registry.tsv")
cat("Registry GWAS:", nrow(reg), "\n")

# Human-readable trait / phenotype family parsed from study_name (data-derived).
trait_of <- function(nm) {
  u <- toupper(nm)
  if (grepl("NASH", u))                         return("NASH/MASH")
  if (grepl("NAFLD|NAFL|MASLD", u))             return("NAFLD/MASLD")
  if (grepl("CIRRHOSIS|CHIRHEP", u))            return("Cirrhosis")
  if (grepl("HCC|HEPATOCELL", u))               return("HCC")
  if (grepl("PDFF", u))                         return("Liver fat (PDFF)")
  if (grepl("ALT", u))                          return("ALT (liver enzyme)")
  if (grepl("AST", u))                          return("AST (liver enzyme)")
  if (grepl("GGT", u))                          return("GGT (liver enzyme)")
  "Other"
}
# Cohort / consortium family parsed from study_name.
cohort_of <- function(nm) {
  if (grepl("FinnGen", nm))            return("FinnGen")
  if (grepl("^BBJ", nm))               return("BioBank Japan")
  if (grepl("PanUKBB", nm))            return("Pan-UKBB")
  if (grepl("^UKBB", nm))              return("UK Biobank")
  if (grepl("deCode|deCODE", nm))      return("deCODE")
  if (grepl("Intermountain", nm))      return("Intermountain")
  if (grepl("Ghouse", nm))             return("Ghouse (UKBB+iPSYCH+HUNT)")
  if (grepl("Sveinbjornsson", nm))     return("Sveinbjornsson meta (deCODE+UKBB+Intermountain+FinnGen)")
  # numeric-prefixed published GWAS carry the PMID after the year
  m <- regmatches(nm, regexpr("^[0-9]{4}_[0-9]+", nm))
  if (length(m) == 1) return(paste0("Published (PMID ", sub("^[0-9]{4}_", "", m), ")"))
  "Other"
}
# SuSiE-coloc reliability is gated on GWAS-side power: <~1000 cases → SuSiE rarely
# converges and coloc.abf is the reported posterior (abf is N-invariant when varbeta
# is supplied, so this is not a loss of information for underpowered cohorts).

reg[, trait        := vapply(study_name, trait_of,  character(1))]
reg[, cohort       := vapply(study_name, cohort_of, character(1))]
reg[, design       := ifelse(trait_type == "binary",
                             sprintf("case-control (%s cases / %s controls)",
                                     format(N_cases, big.mark=","),
                                     format(N_tot - N_cases, big.mark=",")),
                             sprintf("quantitative (N=%s)", format(N_tot, big.mark=",")))]

# ── 2. Fine-mapping coverage (combined_finemapping.csv, per study) ────────────
fm_file <- "results/combined_finemapping.csv"
fm_stats <- data.table(study = character(0))
if (file.exists(fm_file)) {
  fm <- fread(fm_file)
  # coerce logical-ish columns robustly
  to_lgl <- function(x) { if (is.logical(x)) return(x); tolower(as.character(x)) %in% c("true","t","1") }
  fm[, .either_cs := if ("either_in_cs" %in% names(fm)) to_lgl(either_in_cs) else NA]
  fm[, .conv      := if ("susie_converged" %in% names(fm)) to_lgl(susie_converged) else NA]
  pip_col <- if ("recommended_pip" %in% names(fm)) "recommended_pip" else
             if ("max_pip" %in% names(fm)) "max_pip" else NA_character_
  fm_stats <- fm[, .(
    fm_n_variants      = .N,
    fm_n_cs_variants   = sum(.either_cs, na.rm = TRUE),
    fm_n_pip50         = if (!is.na(pip_col)) sum(get(pip_col) > 0.5, na.rm = TRUE) else NA_integer_,
    fm_susie_converged = sum(.conv, na.rm = TRUE)
  ), by = .(study)]
  cat("Fine-mapping studies:", nrow(fm_stats), "\n")
}

# ── 3. COLOC coverage (per-GWAS susie_coloc dir) ─────────────────────────────
coloc_for <- function(study) {
  d <- file.path("results/susie_coloc", study)
  files <- list.files(d, pattern = "^susie_coloc_chr[0-9]+\\.csv$", full.names = TRUE)
  out <- list(coloc_run = FALSE, coloc_n_chr = 0L, coloc_n_genes = 0L,
              coloc_abf_pp4_05 = NA_integer_, coloc_abf_pp4_08 = NA_integer_,
              susie_n_converged = NA_integer_, susie_pp4_05 = NA_integer_,
              susie_pp4_08 = NA_integer_, susie_pp4_09 = NA_integer_,
              coloc_method = "not_run")
  if (length(files) == 0) return(out)
  dt <- rbindlist(lapply(files, function(f) tryCatch(fread(f), error=function(e) NULL)),
                  fill = TRUE)
  if (is.null(dt) || nrow(dt) == 0) return(out)
  out$coloc_run    <- TRUE
  out$coloc_n_chr  <- length(files)
  out$coloc_n_genes <- nrow(dt)
  if ("PP.H4.abf" %in% names(dt)) {
    out$coloc_abf_pp4_05 <- sum(dt[["PP.H4.abf"]] > 0.5, na.rm = TRUE)
    out$coloc_abf_pp4_08 <- sum(dt[["PP.H4.abf"]] > 0.8, na.rm = TRUE)
  }
  if ("method" %in% names(dt)) {
    is_susie <- dt$method == "susie"
    out$susie_n_converged <- sum(is_susie, na.rm = TRUE)
    if ("PP.H4.susie" %in% names(dt)) {
      out$susie_pp4_05 <- sum(is_susie & dt[["PP.H4.susie"]] > 0.5, na.rm = TRUE)
      out$susie_pp4_08 <- sum(is_susie & dt[["PP.H4.susie"]] > 0.8, na.rm = TRUE)
      out$susie_pp4_09 <- sum(is_susie & dt[["PP.H4.susie"]] > 0.9, na.rm = TRUE)
    }
    # classify the analysis actually delivered
    out$coloc_method <- if (out$susie_n_converged > 0) "abf + SuSiE-coloc"
                        else "coloc.abf (SuSiE non-convergent)"
  } else {
    out$coloc_method <- "coloc.abf"
  }
  out
}
coloc_stats <- rbindlist(lapply(reg$study_name, function(s) {
  x <- coloc_for(s); x$study <- s; as.data.table(x)
}), fill = TRUE)

# ── 4. Assemble ───────────────────────────────────────────────────────────────
tab <- merge(reg, fm_stats, by.x = "study_name", by.y = "study", all.x = TRUE)
tab <- merge(tab, coloc_stats, by.x = "study_name", by.y = "study", all.x = TRUE)

# genome build: registry sumstats are all hg19 (reformatted_hg19 / preprocessed = lifted)
tab[, genome_build := "hg19 (GRCh37)"]

# ordering: trait family, then descending N
trait_order <- c("NAFLD/MASLD","NASH/MASH","Cirrhosis","HCC",
                 "Liver fat (PDFF)","ALT (liver enzyme)","AST (liver enzyme)","GGT (liver enzyme)","Other")
tab[, .torder := match(trait, trait_order)]
setorder(tab, .torder, -N_tot)

sel <- tab[, .(
  study_name, cohort, trait, ancestry, design,
  N_tot, N_cases, genome_build, ld_panel,
  fm_n_variants, fm_n_cs_variants, fm_n_pip50, fm_susie_converged,
  coloc_run, coloc_n_chr, coloc_n_genes, coloc_method,
  coloc_abf_pp4_05, coloc_abf_pp4_08,
  susie_n_converged, susie_pp4_05, susie_pp4_08, susie_pp4_09
)]

out_csv <- "results/gwas_analysis_provenance.csv"
fwrite(sel, out_csv)
cat("\nWrote:", out_csv, "(", nrow(sel), "GWAS )\n")

# ── 5. Paper-ready markdown ──────────────────────────────────────────────────
fmt <- function(x) ifelse(is.na(x), "—", format(x, big.mark=",", trim=TRUE))
md <- c(
  "# Supplementary Table SX. GWAS datasets and genetic analyses",
  "",
  sprintf("Generated %s from on-disk artifacts (registry + combined_finemapping.csv + per-GWAS SuSiE-coloc).",
          format(Sys.time(), "%Y-%m-%d %H:%M %Z")),
  "All summary statistics harmonized to hg19/GRCh37. Fine-mapping: SuSiE + CARMA (single-ancestry),",
  "SuSiEx/MESuSiE (multi-ancestry). Colocalization: coloc.abf (primary; N-invariant given varbeta) and",
  "SuSiE-coloc against ancestry-matched LD (PolyFun UKBB-EUR / 1000G EAS/AFR/SAS). PP.H4 = posterior of",
  "shared causal variant. SuSiE-coloc requires GWAS-side fine-mapping convergence, which fails for",
  "underpowered cohorts (<~1,000 cases) — those report coloc.abf only.",
  "Fine-mapping is per-GWAS-locus (SuSiE/CARMA on lead loci); rows marked 'not finemapped' were",
  "colocalized via coloc.abf / multi-ancestry SuSiEx without a single-ancestry SuSiE credible set.",
  "NOTE: regenerate after any in-flight SuSiE-coloc job completes — rows for currently-running GWAS",
  "carry interim/pre-rerun posteriors.",
  "",
  paste("| GWAS | Cohort | Trait | Anc. | Design | N | Finemap vars (CS / PIP>0.5) |",
        "COLOC genes | Method | abf PP4>0.5 (>0.8) | SuSiE conv. | SuSiE PP4>0.5 (>0.8/>0.9) |"),
  paste("|---|---|---|---|---|---|---|---|---|---|---|---|")
)
for (i in seq_len(nrow(sel))) {
  r <- sel[i]
  fm_cell <- if (is.na(r$fm_n_variants)) "not finemapped" else
    sprintf("%s (%s / %s)", fmt(r$fm_n_variants), fmt(r$fm_n_cs_variants), fmt(r$fm_n_pip50))
  coloc_cell <- if (!isTRUE(r$coloc_run)) "—" else fmt(r$coloc_n_genes)
  abf_cell   <- if (is.na(r$coloc_abf_pp4_05)) "—" else
    sprintf("%s (%s)", fmt(r$coloc_abf_pp4_05), fmt(r$coloc_abf_pp4_08))
  sus_conv   <- if (is.na(r$susie_n_converged)) "—" else fmt(r$susie_n_converged)
  sus_cell   <- if (is.na(r$susie_pp4_05) || (!is.na(r$susie_n_converged) && r$susie_n_converged == 0)) "—" else
    sprintf("%s (%s/%s)", fmt(r$susie_pp4_05), fmt(r$susie_pp4_08), fmt(r$susie_pp4_09))
  md <- c(md, sprintf("| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |",
    r$study_name, r$cohort, r$trait, r$ancestry, r$design,
    fmt(r$N_tot), fm_cell, coloc_cell, r$coloc_method, abf_cell, sus_conv, sus_cell))
}
# portfolio totals
md <- c(md, "",
  sprintf("**Portfolio:** %d GWAS across %d ancestries (%s); %d traits. %d GWAS fine-mapped, %d colocalized.",
          nrow(sel), uniqueN(sel$ancestry), paste(sort(unique(sel$ancestry)), collapse="/"),
          uniqueN(sel$trait), sum(!is.na(sel$fm_n_variants)), sum(sel$coloc_run, na.rm=TRUE)),
  sprintf("SuSiE-coloc converged for %d GWAS; coloc.abf-only for %d (underpowered / SuSiE non-convergent).",
          sum(sel$susie_n_converged > 0, na.rm=TRUE),
          sum(sel$coloc_run & (is.na(sel$susie_n_converged) | sel$susie_n_converged == 0), na.rm=TRUE)))
out_md <- "results/gwas_analysis_provenance.md"
writeLines(md, out_md)
cat("Wrote:", out_md, "\n")

# ── 6. Console summary ───────────────────────────────────────────────────────
cat("\n================= PORTFOLIO SUMMARY =================\n")
cat(sprintf("  GWAS total:            %d\n", nrow(sel)))
cat(sprintf("  Ancestries:            %s\n", paste(sort(unique(sel$ancestry)), collapse=", ")))
cat(sprintf("  Trait families:        %s\n", paste(sort(unique(sel$trait)), collapse=", ")))
cat(sprintf("  Fine-mapped:           %d\n", sum(!is.na(sel$fm_n_variants))))
cat(sprintf("  COLOC run:             %d\n", sum(sel$coloc_run, na.rm=TRUE)))
cat(sprintf("  SuSiE-coloc converged: %d\n", sum(sel$susie_n_converged > 0, na.rm=TRUE)))
cat(sprintf("  coloc.abf-only:        %d\n",
            sum(sel$coloc_run & (is.na(sel$susie_n_converged) | sel$susie_n_converged == 0), na.rm=TRUE)))
cat("====================================================\n")
