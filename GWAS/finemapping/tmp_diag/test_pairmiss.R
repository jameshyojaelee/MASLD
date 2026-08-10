# Does the new position-integrity detector reproduce the mode contrast?
# Predicted: `join` (hg19-native) ~0.0005, `join_lift` / `existing` elevated.
suppressPackageStartupMessages(library(data.table))
FM_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
src <- readLines(file.path(FM_DIR, "src/06_susie_coloc.R"))
i0 <- grep("^PANEL_AF_ROOT <- ", src)
closes <- grep("^\\}$", src)
i1 <- min(closes[closes > grep("^strand_verdict_vs_panel <- ", src)])
eval(parse(text = paste(src[i0:i1], collapse = "\n")))
af <- fread(file.path(FM_DIR, "config/gwas_af_sources.tsv"))
CHR <- 22
pa <- load_panel_af("eur", CHR)
tests <- c("2021_34841290_NAFLD_EUR", "2021_34128465_PDFF_EUR",   # join, hg19-native
           "FinnGen_NAFLD", "MVP_ALT_EUR",                        # join_lift
           "2023_36280732_NAFLD_UKBB_EUR")                        # existing
cat(sprintf("%-38s %-10s %10s %10s %s\n","study","mode","n_nonpal","PAIRMISS","guard"))
for (s in tests) {
  h <- af[study_name == s]
  if (!nrow(h) || !file.exists(h$af_sumstats_path[[1]])) next
  g <- fread(h$af_sumstats_path[[1]])[chromosome == CHR]
  if (!nrow(g)) next
  g <- unique(g, by = "position")
  g[pa, on = .(position), `:=`(pg_a1 = i.bim_a1, pg_a2 = i.bim_a2, pg_af = i.af_a1)]
  g[, a1 := toupper(allele1)]; g[, a2 := toupper(allele2)]
  g[, is_pal := (a1 %in% c("A","T") & a2 %in% c("A","T")) |
                (a1 %in% c("C","G") & a2 %in% c("C","G"))]
  np <- g[is_pal == FALSE & is.finite(pg_af)]
  pm <- if (nrow(np) >= 100) mean(is.na(panel_freq_of(np$a1, np$a2, np$pg_a1, np$pg_a2, np$pg_af))) else NA_real_
  cat(sprintf("%-38s %-10s %10d %10.5f %s\n", s, h$mode[[1]], nrow(np), pm,
              if (is.finite(pm) && pm <= 0.02) "ACTIVE" else "DISABLED"))
}
