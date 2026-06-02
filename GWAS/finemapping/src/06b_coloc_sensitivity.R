#!/usr/bin/env Rscript
# 06b_coloc_sensitivity.R
# ---------------------------------------------------------------------------
# p12-prior ROBUSTNESS for the top COLOC hits (standards review: coloc.abf's PP.H4
# depends on the arbitrary p12 prior, default 5e-6; Wallace 2020 recommends reporting
# coloc::sensitivity()). For each target gene's best-GWAS locus we rebuild the coloc.abf
# object exactly as Script 06 does, then recompute PP.H4 across a p12 grid (via the same
# internal coloc:::combine.abf that coloc::sensitivity uses) and report whether the
# colocalization rule (H4 > rule) survives across the field-reasonable p12 window
# [1e-6, 1e-5] around the 5e-6 default. Robust hits survive; fragile hits flip with p12.
#
# Targets (default): genes whose coloc_best_pp4 OR coloc_best_susie_pp4 > THRESH in
# gene_level_coloc.csv. Override with COLOC_SENS_GENES="GENE1,GENE2,...".
# Output: results/susie_coloc/coloc_sensitivity.csv  (+ console summary)
# Run AFTER the COLOC re-run + 07 aggregate. Light: micromamba run -n rnaseq Rscript src/06b_coloc_sensitivity.R
# ---------------------------------------------------------------------------
suppressMessages({ library(data.table); library(coloc) })

BASE_DIR <- Sys.getenv("MASLD_PROJECT_ROOT",
  unset = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FM_DIR   <- file.path(BASE_DIR, "GWAS/finemapping")
setwd(FM_DIR)
EQTL_DIR <- file.path(BASE_DIR, "data/broadaway_eqtl")
COLOC_DIR <- file.path(FM_DIR, "results/susie_coloc")

# Match Script 06 exactly
COLOC_P1 <- 1e-4; COLOC_P2 <- 1e-4; COLOC_P12 <- 5e-6; EQTL_N <- 1183; MIN_SNPS <- 100L
THRESH    <- as.numeric(Sys.getenv("COLOC_SENS_THRESH", "0.5"))   # which genes count as "top hits"
RULE      <- as.numeric(Sys.getenv("COLOC_SENS_RULE",   "0.5"))   # H4 > RULE is the coloc claim
P12_LO <- 1e-6; P12_HI <- 1e-5                                     # field-reasonable window around 5e-6
P12_GRID <- 10^seq(-7, -4, length.out = 31)

registry <- read.delim("config/gwas_registry.tsv", stringsAsFactors = FALSE)

# ---- target gene list ----------------------------------------------------
gl <- fread(file.path(COLOC_DIR, "gene_level_coloc.csv"))
gcol <- intersect(c("gene", "gene_name", "gene_symbol"), names(gl))[1]
pp_abf <- if ("coloc_best_pp4" %in% names(gl)) gl$coloc_best_pp4 else rep(NA, nrow(gl))
pp_sus <- if ("coloc_best_susie_pp4" %in% names(gl)) gl$coloc_best_susie_pp4 else rep(NA, nrow(gl))
bestg  <- if ("coloc_best_gwas" %in% names(gl)) gl$coloc_best_gwas else
          if ("coloc_best_susie_gwas" %in% names(gl)) gl$coloc_best_susie_gwas else NA
gl[, `:=`(.sym = get(gcol), .pp = pmax(pp_abf, pp_sus, na.rm = TRUE), .bg = bestg)]

env_genes <- Sys.getenv("COLOC_SENS_GENES", "")
if (nzchar(env_genes)) {
  want <- trimws(strsplit(env_genes, ",")[[1]])
  targets <- gl[.sym %in% want]
} else {
  targets <- gl[.pp > THRESH]
}
targets <- targets[!is.na(.bg) & .bg != "" & .bg != "NA"]
targets[, .bg := gsub('"', '', .bg)]
cat(sprintf("Targets: %d genes (PP>%.2f or explicit list), across %d GWAS\n",
            nrow(targets), THRESH, length(unique(targets$.bg))))
if (!nrow(targets)) { cat("No targets — nothing to do.\n"); quit(status = 0) }

# ---- symbol -> (ENSG, chr) map from the eQTL files (cheap: 3 cols) --------
cat("Building gene->chr map from Broadaway eQTL headers...\n")
sym_map <- rbindlist(lapply(1:22, function(ch) {
  f <- file.path(EQTL_DIR, sprintf("chr%d_marginal_summary_results.tsv", ch))
  if (!file.exists(f)) return(NULL)
  d <- unique(fread(f, select = c("ENSG", "GeneSymbol"))); d[, chr := ch]; d
}), fill = TRUE)
targets <- merge(targets, sym_map, by.x = ".sym", by.y = "GeneSymbol", all.x = TRUE)
targets <- targets[!is.na(chr)]
cat(sprintf("  mapped %d targets to a chromosome\n", nrow(targets)))

# ---- harmonize one gene's GWAS x eQTL locus (mirrors Script 06) -----------
build_merged <- function(eqtl_gene, gwas) {
  eqtl_gene <- copy(eqtl_gene); eqtl_gene[, merge_key := paste(CHR, POS, sep = ":")]
  m <- merge(eqtl_gene[, .(merge_key, eqtl_ea = EA, eqtl_nea = NEA,
                           eqtl_beta = Beta, eqtl_se = SE)],
             gwas[, .(merge_key, gwas_a1 = allele1, gwas_a2 = allele2,
                      gwas_beta = beta, gwas_se = se, gwas_pval = pval)], by = "merge_key")
  m <- m[order(gwas_pval)][!duplicated(merge_key)]
  if (!nrow(m)) return(m)
  m[, am := (eqtl_ea == gwas_a1 & eqtl_nea == gwas_a2)]
  m[, af := (eqtl_ea == gwas_a2 & eqtl_nea == gwas_a1)]
  m <- m[am | af]; m[af == TRUE, eqtl_beta := -eqtl_beta]
  m[!((gwas_a1 %in% c("A","T") & gwas_a2 %in% c("A","T")) |
      (gwas_a1 %in% c("C","G") & gwas_a2 %in% c("C","G")))]
}

pp4_over_p12 <- function(abf_res, p12) {
  cb <- NULL
  invisible(capture.output(suppressMessages(   # combine.abf prints via cat(); silence it
    cb <- coloc:::combine.abf(abf_res$results$lABF.df1, abf_res$results$lABF.df2,
                              p1 = COLOC_P1, p2 = COLOC_P2, p12 = p12))))
  as.numeric(cb["PP.H4.abf"])
}

# ---- loop by (GWAS, chr): load once, process its target genes -------------
out <- list()
for (g in unique(targets$.bg)) {
  row <- registry[registry$study_name == g, ]
  if (!nrow(row)) { cat("  [skip] GWAS not in registry:", g, "\n"); next }
  if (!file.exists(row$sumstats_path)) { cat("  [skip] sumstats missing:", g, "\n"); next }
  gtype <- ifelse(row$trait_type == "binary", "cc", "quant")
  for (ch in sort(unique(targets[.bg == g, chr]))) {
    eqf <- file.path(EQTL_DIR, sprintf("chr%d_marginal_summary_results.tsv", ch))
    if (!file.exists(eqf)) next
    gw <- fread(row$sumstats_path); gw <- gw[chromosome == ch]
    gw[, merge_key := paste(chromosome, position, sep = ":")]
    eq <- fread(eqf)
    for (gene in targets[.bg == g & chr == ch, .sym]) {
      ensg <- targets[.bg == g & chr == ch & .sym == gene, ENSG][1]
      eg <- eq[ENSG == ensg]
      if (nrow(eg) < MIN_SNPS) next
      m <- build_merged(eg, gw)
      if (nrow(m) < MIN_SNPS) next
      res <- tryCatch({
        d1 <- list(beta = m$gwas_beta, varbeta = m$gwas_se^2, N = row$N_tot,
                   type = gtype, snp = m$merge_key)
        if (gtype == "cc" && row$N_cases > 0) d1$s <- row$N_cases / row$N_tot
        if (gtype == "quant") d1$sdY <- 1
        d2 <- list(beta = m$eqtl_beta, varbeta = m$eqtl_se^2, N = EQTL_N,
                   type = "quant", sdY = 1, snp = m$merge_key)
        suppressMessages(suppressWarnings(coloc.abf(d1, d2, p1 = COLOC_P1, p2 = COLOC_P2, p12 = COLOC_P12)))
      }, error = function(e) NULL)
      if (is.null(res)) next
      grid_pp4 <- vapply(P12_GRID, function(p) pp4_over_p12(res, p), numeric(1))
      pp_def <- pp4_over_p12(res, COLOC_P12)
      pp_lo  <- pp4_over_p12(res, P12_LO)
      pp_hi  <- pp4_over_p12(res, P12_HI)
      # smallest/largest p12 in grid where rule holds
      passing <- P12_GRID[grid_pp4 > RULE]
      out[[length(out) + 1L]] <- data.table(
        gene = gene, ensembl = ensg, chr = ch, gwas = g, n_snps = nrow(m),
        pp4_default = round(pp_def, 4),
        pp4_at_1e6 = round(pp_lo, 4), pp4_at_1e5 = round(pp_hi, 4),
        p12_pass_min = if (length(passing)) min(passing) else NA_real_,
        p12_pass_max = if (length(passing)) max(passing) else NA_real_,
        robust = (pp_lo > RULE & pp_hi > RULE),                 # holds across [1e-6,1e-5]
        rule = paste0("H4>", RULE))
    }
  }
}
res <- rbindlist(out, fill = TRUE)
if (!nrow(res)) { cat("No coloc objects rebuilt (targets may lack overlap). Done.\n"); quit(status = 0) }
setorder(res, -pp4_default)
of <- file.path(COLOC_DIR, "coloc_sensitivity.csv")
fwrite(res, of)
cat(sprintf("\nWritten: %s (%d genes)\n", of, nrow(res)))
cat(sprintf("ROBUST (H4>%.1f across p12 in [1e-6,1e-5]): %d / %d\n",
            RULE, sum(res$robust, na.rm = TRUE), nrow(res)))
cat("\n=== fragile hits (flip within the p12 window) ===\n")
print(res[robust == FALSE, .(gene, gwas, pp4_at_1e6, pp4_default, pp4_at_1e5)])
cat("\n=== top robust hits ===\n")
print(head(res[robust == TRUE, .(gene, gwas, pp4_default, pp4_at_1e6, pp4_at_1e5)], 12))
