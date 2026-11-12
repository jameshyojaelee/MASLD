#!/usr/bin/env Rscript
# ENZFAT1 S2 step 2: per-locus liver-fat and cirrhosis effect at the 169 ALT loci.
# Unit: loci. All effects are aligned to the ALT-raising allele in Pan-UKBB EUR ALT.
# ashr (adaptive shrinkage; built for a set of effect estimates with standard errors) per outcome.
# Family F3: lfsr < 0.05 over the loci for liver fat; cirrhosis is secondary.
# Usage: 21_locus_effects.R <locus_rows.tsv> <loci.tsv> <out dir>
suppressMessages({library(data.table); library(ashr)})
args <- commandArgs(trailingOnly = TRUE)
rows <- fread(args[1]); loci <- fread(args[2]); out <- args[3]
set.seed(20261110)
rows[, locus := as.character(locus)]; loci[, locus := as.character(locus)]

ref <- rows[trait == "alt_ukb", .(locus, ref_allele = ifelse(beta > 0, ea, oa), alt_beta = abs(beta), alt_se = se)]
d <- merge(rows, ref, by = "locus")
d[, beta_al := ifelse(ea == ref_allele, beta, -beta)]
stopifnot(d[trait == "alt_ukb", all(beta_al > 0)])
wide <- dcast(d, locus + ref_allele ~ trait, value.var = c("beta_al", "se"))
info <- loci[, .(locus, chr, pos, genes, cls_DIG1 = cls)]
wide <- merge(info, wide, by = "locus", all.y = TRUE)
is_locus <- !grepl("^ctrl_", wide$locus)

# adaptive shrinkage per outcome over the lead variants of the locus table (controls excluded from the fit)
for (tr in c("pdff", "cirr_ghouse", "nafld_finngen")) {
  b <- wide[[paste0("beta_al_", tr)]]; s <- wide[[paste0("se_", tr)]]
  use <- is_locus & is.finite(b) & is.finite(s)
  fit <- ash(b[use], s[use], mixcompdist = "halfuniform")
  wide[use, paste0(c("postmean_", "postsd_", "lfsr_", "lfdr_"), tr) :=
         .(get_pm(fit), get_psd(fit), get_lfsr(fit), get_lfdr(fit))]
  # control rows outside the locus table: posterior under the prior fitted on the loci (g fixed)
  uc <- !is_locus & is.finite(b) & is.finite(s)
  if (any(uc)) {
    fc <- ash(b[uc], s[uc], g = get_fitted_g(fit), fixg = TRUE, mixcompdist = "halfuniform")
    wide[uc, paste0(c("postmean_", "postsd_", "lfsr_", "lfdr_"), tr) :=
           .(get_pm(fc), get_psd(fc), get_lfsr(fc), get_lfdr(fc))]
  }
  cat(tr, ": loci fitted", sum(use), "; pi0", round(get_pi0(fit), 3), "; lfsr < 0.05:",
      sum(get_lfsr(fit) < 0.05), "\n")
}
# effect per unit ALT effect, delta-method SE (ALT estimate from Pan-UKBB; sensitivity with MVP)
ratio <- function(f, sf, a, sa) list(f / a, sqrt(sf^2 / a^2 + f^2 * sa^2 / a^4))
wide[, c("fat_per_alt", "fat_per_alt_se") := ratio(beta_al_pdff, se_pdff, beta_al_alt_ukb, se_alt_ukb)]
wide[, c("cirr_per_alt", "cirr_per_alt_se") := ratio(beta_al_cirr_ghouse, se_cirr_ghouse, beta_al_alt_ukb, se_alt_ukb)]
wide[, c("fat_per_alt_mvp", "fat_per_alt_mvp_se") := ratio(beta_al_pdff, se_pdff, beta_al_alt_mvp, se_alt_mvp)]
for (tr in c("pdff", "cirr_ghouse", "nafld_finngen", "nafld_mvp", "cirr_mvp", "albumin_mvp", "platelet_mvp", "alt_mvp"))
  wide[, paste0("z_", tr) := get(paste0("beta_al_", tr)) / get(paste0("se_", tr))]
fwrite(wide, file.path(out, "locus_effects.tsv"), sep = "\t")

# controls
has <- function(g) sapply(strsplit(wide$genes, "[,;]"), function(x) g %in% x)
gpt <- wide[chr == 8 & abs(pos - 145733794) < 5e5]
ctl <- rbind(
  # PNPLA3 is a control row (rs738409); TM6SF2 is the locus whose lead is rs58542926 (19:19379549);
  # the locus table lists neither gene name, so they are found by variant, not by gene
  wide[locus == "ctrl_PNPLA3_rs738409" | (chr == 19 & pos == 19379549) | has("GCKR"), .(control = "liver-fat effect expected (lfsr < 0.05)", locus, genes, z_pdff, lfsr_pdff, z_cirr_ghouse, pass = lfsr_pdff < 0.05)],
  gpt[, .(control = "GPT locus: liver-fat |z| < 2", locus, genes, z_pdff, lfsr_pdff, z_cirr_ghouse, pass = abs(z_pdff) < 2)],
  wide[has("GOT1") | locus == "ctrl_GOT1_lead", .(control = "GOT1 locus: liver-fat |z| < 2", locus, genes, z_pdff, lfsr_pdff, z_cirr_ghouse, pass = abs(z_pdff) < 2)],
  wide[grepl("^ctrl_HSD17B13", locus), .(control = "HSD17B13: ALT-raising allele raises cirrhosis, liver-fat |z| < 2", locus, genes, z_pdff, lfsr_pdff, z_cirr_ghouse, pass = abs(z_pdff) < 2 & z_cirr_ghouse > 0)])
fwrite(ctl, file.path(out, "controls.tsv"), sep = "\t")
print(ctl)
# ascertainment check at the GPT locus: MVP 571.5 z over MVP ALT z, beside the LDSC cross-trait intercept
if (nrow(gpt)) cat("GPT locus z(MVP 571.5) / z(MVP ALT):", round(gpt$z_nafld_mvp / gpt$z_alt_mvp, 3), "\n")
# the GPT lead is an indel absent from the liver-fat file, so it has no rsID for the MVP files; report the
# same ratio from the MVP z scores already in the DIG1 locus table (aligned to the ALT-raising allele there)
g68 <- loci[chr == 8 & abs(pos - 145733794) < 5e5]
if (nrow(g68)) cat("GPT locus, DIG1 table z(MVP 571.5) / z(MVP ALT):", round(g68$z_MVP_NAFLD_EUR / g68$z_MVP_ALT_EUR, 3), "\n")
summ <- data.table(
  quantity = c("loci with liver-fat and ALT rows", "liver fat lfsr < 0.05 (F3)", "cirrhosis lfsr < 0.05 (secondary)",
               "FinnGen NAFLD lfsr < 0.05 (out of sample)", "both liver fat and cirrhosis lfsr < 0.05",
               "median fat effect per unit ALT effect"),
  value = c(sum(is_locus & is.finite(wide$lfsr_pdff)), sum(is_locus & wide$lfsr_pdff < 0.05, na.rm = TRUE),
            sum(is_locus & wide$lfsr_cirr_ghouse < 0.05, na.rm = TRUE), sum(is_locus & wide$lfsr_nafld_finngen < 0.05, na.rm = TRUE),
            sum(is_locus & wide$lfsr_pdff < 0.05 & wide$lfsr_cirr_ghouse < 0.05, na.rm = TRUE),
            median(wide$fat_per_alt[is_locus], na.rm = TRUE)))
fwrite(summ, file.path(out, "summary.tsv"), sep = "\t"); print(summ)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
cat("S2 DONE\n")
