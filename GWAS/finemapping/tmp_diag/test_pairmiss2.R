# Does the WINDOWED detector reproduce the mode contrast the chromosome-wide
# one missed? Predicted: opposite-calls sit in high-pairmiss windows for
# join_lift/existing, and NOT for hg19-native join.
suppressPackageStartupMessages(library(data.table))
FM_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
src <- readLines(file.path(FM_DIR, "src/06_susie_coloc.R"))
i0 <- grep("^PANEL_AF_ROOT <- ", src)
closes <- grep("^\\}$", src)
i1 <- min(closes[closes > grep("^strand_verdict_vs_panel <- ", src)])
eval(parse(text = paste(src[i0:i1], collapse = "\n")))
af <- fread(file.path(FM_DIR, "config/gwas_af_sources.tsv")); CHR <- 22
pa <- load_panel_af("eur", CHR)
PM_MAX <- 0.02; PM_MIN_N <- 20L; W <- 1e5
cat(sprintf("%-38s %-10s %8s %10s %10s %s\n","study","mode","n_opp","medPM_opp","medPM_same","opp windows blocked"))
for (s in c("2021_34841290_NAFLD_EUR","2021_34128465_PDFF_EUR","FinnGen_NAFLD",
            "MVP_ALT_EUR","2023_36280732_NAFLD_UKBB_EUR")) {
  h <- af[study_name == s]; if (!nrow(h) || !file.exists(h$af_sumstats_path[[1]])) next
  g <- unique(fread(h$af_sumstats_path[[1]])[chromosome == CHR], by = "position")
  if (!nrow(g)) next
  g[pa, on = .(position), `:=`(pg_a1 = i.bim_a1, pg_a2 = i.bim_a2, pg_af = i.af_a1)]
  g[, a1 := toupper(allele1)][, a2 := toupper(allele2)]
  g[, is_pal := (a1 %in% c("A","T") & a2 %in% c("A","T")) | (a1 %in% c("C","G") & a2 %in% c("C","G"))]
  g[, np_miss := ifelse(is_pal == FALSE & is.finite(pg_af),
        as.numeric(is.na(panel_freq_of(a1, a2, pg_a1, pg_a2, pg_af))), NA_real_)]
  g[, win := floor(position / W)]
  g[, win_pm := { n <- sum(!is.na(np_miss)); if (n >= PM_MIN_N) rep(mean(np_miss, na.rm=TRUE), .N) else rep(NA_real_, .N) }, by = win]
  g[, pf := panel_freq_of(a1, a2, pg_a1, pg_a2, pg_af)]
  g[, verdict := strand_verdict_vs_panel(af, pf)]
  o <- g[is_pal == TRUE & verdict == "opposite"]; sm <- g[is_pal == TRUE & verdict == "same"]
  cat(sprintf("%-38s %-10s %8d %10.5f %10.5f %d of %d\n", s, h$mode[[1]], nrow(o),
      median(o$win_pm, na.rm=TRUE), median(sm$win_pm, na.rm=TRUE),
      sum(!(is.finite(o$win_pm) & o$win_pm <= PM_MAX)), nrow(o)))
}
