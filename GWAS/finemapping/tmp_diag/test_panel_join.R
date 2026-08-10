#!/usr/bin/env Rscript
# Exercise the EXACT helpers now in 06_susie_coloc.R against real chr22 data.
# Purpose: prove the join finds variants (not silently all-NA) and that the
# verdicts behave -- including PNPLA3 rs738409, the motivating case.
suppressPackageStartupMessages(library(data.table))
FM_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
src <- readLines(file.path(FM_DIR, "src/06_susie_coloc.R"))
i0 <- grep("^PANEL_AF_ROOT <- ", src); i1 <- grep("^\\}$", src); i1 <- min(i1[i1 > grep("^strand_verdict_vs_panel <- ", src)])
eval(parse(text = paste(src[i0:i1], collapse = "\n")))
cat("helpers loaded from 06_susie_coloc.R lines", i0, "-", i1, "\n\n")

pa_eur <- load_panel_af("eur", 22); pa_afr <- load_panel_af("afr", 22)
cat(sprintf("panel eur chr22: %d variants | afr chr22: %d\n\n", nrow(pa_eur), nrow(pa_afr)))

# --- PNPLA3 rs738409, chr22:44324727 (hg19), C/G palindrome ---
r <- pa_eur[position == 44324727]
cat("=== PNPLA3 rs738409 ===\n"); print(r)
eqtl_ea <- "C"; eqtl_eaf <- 0.7732   # Broadaway, per upstream trace
pf <- panel_freq_of(eqtl_ea, r$bim_a1, r$bim_a2, r$af_a1)
cat(sprintf("  panel freq of eQTL EA (%s) = %.4f   eQTL EAF = %.4f\n", eqtl_ea, pf, eqtl_eaf))
cat(sprintf("  verdict = %s   <- must be 'same' (variant RETAINED)\n\n",
            strand_verdict_vs_panel(eqtl_eaf, pf)))

# --- behavioural check: a deliberately flipped copy must read 'opposite' ---
set.seed(3); s <- pa_eur[is.finite(af_a1)][sample(.N, 20000)]
obs  <- s$af_a1                     # truthful: same strand
flip <- 1 - s$af_a1                 # simulated reverse strand
v_ok <- strand_verdict_vs_panel(obs,  s$af_a1)
v_fl <- strand_verdict_vs_panel(flip, s$af_a1)
cat("=== behaviour on 20,000 real chr22 variants ===\n")
cat(sprintf("  truthful copy : same=%6d  opposite=%6d  abstain=%6d\n",
            sum(v_ok %in% "same"), sum(v_ok %in% "opposite"), sum(is.na(v_ok))))
cat(sprintf("  flipped copy  : same=%6d  opposite=%6d  abstain=%6d\n",
            sum(v_fl %in% "same"), sum(v_fl %in% "opposite"), sum(is.na(v_fl))))
cat(sprintf("  FALSE 'opposite' on truthful data: %d  (must be 0)\n", sum(v_ok %in% "opposite")))
cat(sprintf("  recall on flipped data: %.4f\n", mean(v_fl %in% "opposite")))
