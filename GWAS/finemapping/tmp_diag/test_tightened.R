suppressPackageStartupMessages(library(data.table))
FM_DIR <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
FM <- "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
src <- readLines(file.path(FM, "src/06_susie_coloc.R"))
i0 <- grep("^PANEL_AF_ROOT <- ", src)
closes <- grep("^\\}$", src)
i1 <- min(closes[closes > grep("^strand_verdict_vs_panel <- ", src)])
eval(parse(text = paste(src[i0:i1], collapse = "\n")))
cat("helpers sourced from lines", i0, "-", i1, "\n")
pa <- load_panel_af("eur", 22)
r  <- pa[position == 44324727]
cat(sprintf("PNPLA3 rs738409: panel %s/%s af_a1=%.4f\n", r$bim_a1, r$bim_a2, r$af_a1))
f_ok  <- panel_freq_of("C", "G", r$bim_a1, r$bim_a2, r$af_a1)   # correct pair
f_bad <- panel_freq_of("C", "T", r$bim_a1, r$bim_a2, r$af_a1)   # wrong partner
cat(sprintf("  correct pair C/G : %.4f  verdict=%s  <- must be 'same' (retained)\n",
            f_ok, strand_verdict_vs_panel(0.7732, f_ok)))
cat(sprintf("  wrong   pair C/T : %s  <- must be NA (pair guard fires)\n", f_bad))
stopifnot(is.na(f_bad), !is.na(f_ok))
set.seed(3); s <- pa[is.finite(af_a1)][sample(.N, 20000)]
pf   <- panel_freq_of(s$bim_a1, s$bim_a2, s$bim_a1, s$bim_a2, s$af_a1)
v_ok <- strand_verdict_vs_panel(s$af_a1,     pf)
v_fl <- strand_verdict_vs_panel(1 - s$af_a1, pf)
cat(sprintf("20k variants: false-'opposite' on truthful = %d (must be 0) | recall on flipped = %.4f\n",
            sum(v_ok %in% "opposite"), mean(v_fl %in% "opposite")))
stopifnot(sum(v_ok %in% "opposite") == 0)
cat("TIGHTENED HELPERS OK\n")
