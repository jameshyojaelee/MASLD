#!/usr/bin/env Rscript
# DIR6 step 2 (PRESPEC_DIR6.md sections 4 to 8): per-gene classes, counts, controls, extra columns.
# Unit of inference: independent loci (1 Mb single linkage of instrument positions). Genes reported beside.
# Usage: 02_classify.R <project_root> <step1 out dir>
suppressMessages({library(data.table); library(ashr)})
args <- commandArgs(trailingOnly = TRUE)
ROOT <- args[1]; OUT <- args[2]
set.seed(20261111)
rows <- fread(file.path(OUT, "rows.tsv"), colClasses = list(character = c("id", "trait", "ea", "oa")))
inst <- fread(file.path(OUT, "instruments.tsv"), colClasses = list(character = c("chr", "vid", "ensembl")))
ctrlv <- fread(file.path(OUT, "controls_variants.tsv"), colClasses = list(character = "chr"))
h26 <- fread(file.path(ROOT, "Analysis/Hypotheses/H26-20261007T210504Z/run_22092011/h26_gene_direction_table.tsv"))
enz <- fread(file.path(ROOT, "Analysis/Hypotheses/ENZFAT1-20261110T204401Z/out_S2/locus_effects.tsv"))
stopifnot(nrow(inst) == 286)

# ---------------------------------------------------------------- alignment
# instruments: the H26 trait-raising allele (third field of vid); controls: Pan-UKBB ALT-raising allele
# (GOT1: AST-raising allele)
al <- unique(inst[, .(id = vid, allele = tstrsplit(vid, ":")[[3]])])
for (cid in ctrlv$id) {
  tr <- if (grepl("GOT1", cid)) "AST" else "ALT"
  r <- rows[id == cid & trait == tr]
  if (nrow(r) == 1) al <- rbind(al, data.table(id = cid, allele = if (r$beta > 0) r$ea else r$oa))
}
d <- merge(rows, al, by = "id")
d[, beta_al := fifelse(ea == allele, beta, -beta)]
w <- dcast(d, id ~ trait, value.var = c("beta_al", "se"))
src <- dcast(d[trait %in% c("pdff", "cirr")], id ~ trait, value.var = c("source", "proxy_r2"))
w <- merge(w, src, by = "id", all.x = TRUE)
is_ctrl <- grepl("^ctrl_", w$id)

# ---------------------------------------------------------------- ashr per outcome (unique instrument rows)
fits <- list()
for (tr in c("pdff", "cirr")) {
  b <- w[[paste0("beta_al_", tr)]]; s <- w[[paste0("se_", tr)]]
  use <- !is_ctrl & is.finite(b) & is.finite(s)
  fit <- ash(b[use], s[use], mixcompdist = "halfuniform")
  fits[[tr]] <- fit
  w[use, paste0(c("pm_", "psd_", "lfsr_"), tr) := .(get_pm(fit), get_psd(fit), get_lfsr(fit))]
  uc <- is_ctrl & is.finite(b) & is.finite(s)
  if (any(uc)) {
    fc <- ash(b[uc], s[uc], g = get_fitted_g(fit), fixg = TRUE, mixcompdist = "halfuniform")
    w[uc, paste0(c("pm_", "psd_", "lfsr_"), tr) := .(get_pm(fc), get_psd(fc), get_lfsr(fc))]
  }
  cat(tr, ": instrument rows fitted", sum(use), "; pi0", round(get_pi0(fit), 3), "; lfsr < 0.05:", sum(get_lfsr(fit) < 0.05), "\n")
  w[, paste0("call_", tr) := fifelse(!is.finite(get(paste0("lfsr_", tr))), "untestable",
                              fifelse(get(paste0("lfsr_", tr)) >= 0.05, "no outcome support",
                              fifelse(get(paste0("pm_", tr)) > 0, "concordant", "reversed")))]
}
fwrite(w, file.path(OUT, "instrument_effects.tsv"), sep = "\t")

# ---------------------------------------------------------------- gene table
g <- merge(inst, w[!is_ctrl], by.x = "vid", by.y = "id", all.x = TRUE)
for (tr in c("pdff", "cirr")) g[is.na(get(paste0("call_", tr))), paste0("call_", tr) := "untestable"]
g[, combined := fcase(
  (call_pdff == "concordant" | call_cirr == "concordant") & !(call_pdff == "reversed" | call_cirr == "reversed"), "concordant",
  (call_pdff == "reversed" | call_cirr == "reversed") & !(call_pdff == "concordant" | call_cirr == "concordant"), "reversed",
  (call_pdff == "concordant" | call_cirr == "concordant"), "mixed",
  call_pdff == "untestable" & call_cirr == "untestable", "untestable",
  default = "enzyme support only")]
g[, one_outcome_untestable := combined == "enzyme support only" & (call_pdff == "untestable" | call_cirr == "untestable")]

# enzyme denominator: Pan-UKBB EUR effect of the defining enzyme (ALT for NAFLD/PDFF-defined genes)
g[, denom_trait := fifelse(defining_trait %in% c("ALT", "AST", "GGT"), defining_trait, "ALT")]
g[, enz_beta := fcase(denom_trait == "ALT", beta_al_ALT, denom_trait == "AST", beta_al_AST, denom_trait == "GGT", beta_al_GGT)]
g[, enz_se := fcase(denom_trait == "ALT", se_ALT, denom_trait == "AST", se_AST, denom_trait == "GGT", se_GGT)]
g[, enzyme_sign_differs_in_panukbb_eur := is.finite(enz_beta) & enz_beta <= 0]
g[, enzyme_missing_in_panukbb := !is.finite(enz_beta)]

# equivalence bounds (section 4)
lc <- enz[!grepl("^ctrl_", locus)]
fat_med <- median(lc$fat_per_alt, na.rm = TRUE); cir_med <- median(lc$cirr_per_alt, na.rm = TRUE)
bound_rule_cirr <- "0.5 x median signed ratio"
if (!is.finite(cir_med) || cir_med <= 0) { cir_med <- median(abs(lc$cirr_per_alt), na.rm = TRUE); bound_rule_cirr <- "0.5 x median absolute ratio (signed median <= 0)" }
bounds <- c(pdff = 0.5 * fat_med, cirr = 0.5 * cir_med)
cat("equivalence bounds: liver fat", round(bounds["pdff"], 4), "; cirrhosis", round(bounds["cirr"], 4), "(", bound_rule_cirr, ")\n")

ratio_ci <- function(f, sf, a, sa) { r <- f / a; s <- sqrt(sf^2 / a^2 + f^2 * sa^2 / a^4); list(r, r - 1.959964 * s, r + 1.959964 * s) }
for (tr in c("pdff", "cirr")) {
  g[, paste0(c("ratio_", "ratio_lo_", "ratio_hi_"), tr) := ratio_ci(get(paste0("beta_al_", tr)), get(paste0("se_", tr)), enz_beta, enz_se)]
  g[, paste0("power_", tr) := fifelse(get(paste0("call_", tr)) != "no outcome support", NA_character_,
        fifelse(enzyme_sign_differs_in_panukbb_eur | enzyme_missing_in_panukbb, "not interpretable",
        fifelse(get(paste0("ratio_lo_", tr)) > -bounds[tr] & get(paste0("ratio_hi_", tr)) < bounds[tr],
                "evidence of no effect larger than the bound", "underpowered (no evidence)")))]
}

# independent loci: 1 Mb single linkage on GRCh37 instrument positions
pos <- unique(g[, .(chr, pos37)])[order(as.integer(chr), pos37)]
pos[, locus := cumsum(c(TRUE, chr[-1] != chr[-.N] | diff(pos37) > 1e6))]
g <- merge(g, pos, by = c("chr", "pos37"))
g[, locus_id := paste0("L", locus)]
cat("independent loci:", uniqueN(g$locus_id), "\n")

# ---------------------------------------------------------------- counts
classes <- c("concordant", "reversed", "mixed", "enzyme support only", "untestable")
ci <- function(x, n) if (n > 0) binom.test(x, n)$conf.int else c(NA, NA)
count_tab <- function(dd, label) {
  out <- list()
  for (cl in c("all", "inhibit", "activate")) {
    s <- if (cl == "all") dd else dd[call == cl]
    # genes
    n <- nrow(s)
    for (k in c(classes)) { x <- sum(s$combined == k); out[[length(out) + 1]] <- data.table(set = label, call = cl, unit = "gene", class = k, count = x, n = n, prop = x / n, lo95 = ci(x, n)[1], hi95 = ci(x, n)[2]) }
    for (k in c("concordant on liver fat", "concordant on cirrhosis", "reversed on liver fat", "reversed on cirrhosis", "untestable on liver fat", "untestable on cirrhosis")) {
      x <- switch(k, "concordant on liver fat" = sum(s$call_pdff == "concordant"), "concordant on cirrhosis" = sum(s$call_cirr == "concordant"),
                  "reversed on liver fat" = sum(s$call_pdff == "reversed"), "reversed on cirrhosis" = sum(s$call_cirr == "reversed"),
                  "untestable on liver fat" = sum(s$call_pdff == "untestable"), "untestable on cirrhosis" = sum(s$call_cirr == "untestable"))
      out[[length(out) + 1]] <- data.table(set = label, call = cl, unit = "gene", class = k, count = x, n = n, prop = x / n, lo95 = ci(x, n)[1], hi95 = ci(x, n)[2])
    }
    # loci: class among this call's genes at the locus; heterogeneous if they differ
    L <- s[, .(class = if (uniqueN(combined) == 1) combined[1] else "heterogeneous locus",
               cf = any(call_pdff == "concordant"), cc = any(call_cirr == "concordant"),
               rf = any(call_pdff == "reversed"), rc = any(call_cirr == "reversed")), by = locus_id]
    n <- nrow(L)
    for (k in c(classes, "heterogeneous locus")) { x <- sum(L$class == k); out[[length(out) + 1]] <- data.table(set = label, call = cl, unit = "independent locus", class = k, count = x, n = n, prop = x / n, lo95 = ci(x, n)[1], hi95 = ci(x, n)[2]) }
    for (k in c("any gene concordant on liver fat", "any gene concordant on cirrhosis", "any gene reversed on liver fat", "any gene reversed on cirrhosis")) {
      x <- switch(k, "any gene concordant on liver fat" = sum(L$cf), "any gene concordant on cirrhosis" = sum(L$cc),
                  "any gene reversed on liver fat" = sum(L$rf), "any gene reversed on cirrhosis" = sum(L$rc))
      out[[length(out) + 1]] <- data.table(set = label, call = cl, unit = "independent locus", class = k, count = x, n = n, prop = x / n, lo95 = ci(x, n)[1], hi95 = ci(x, n)[2])
    }
  }
  rbindlist(out)
}
counts <- rbind(count_tab(g, "all 286 calls"),
                count_tab(g[defining_trait %in% c("ALT", "AST", "GGT")], "sensitivity: enzyme-defined only (NAFLD/PDFF-defined removed)"),
                count_tab(g[primary_evidence_class == "genetically_anchored"], "sensitivity: 269 genetically anchored only"))
fwrite(counts, file.path(OUT, "class_counts.tsv"), sep = "\t")
pw <- rbind(g[combined == "enzyme support only", .N, by = .(power = power_pdff)][, outcome := "liver fat"],
            g[combined == "enzyme support only", .N, by = .(power = power_cirr)][, outcome := "cirrhosis"])
pw[, bound := fifelse(outcome == "liver fat", bounds["pdff"], bounds["cirr"])]
fwrite(pw, file.path(OUT, "power_enzyme_only.tsv"), sep = "\t")

# ---------------------------------------------------------------- co-effects of the enzyme-lowering allele
co <- d[trait %in% c("tg_glgc", "bmi_yengo", "homair_manning") & !grepl("^ctrl_", id)]
co[, `:=`(beta_lower = -beta_al, z = -beta_al / se)]
co[, p := 2 * pnorm(-abs(z))]
co[, q := p.adjust(p, "BH")]   # family: all finite variant x trait co-effect tests (TG, BMI, HOMA-IR)
co[, flag := fifelse(q < 0.05, fifelse(beta_lower > 0, "enzyme-lowering allele raises (q < 0.05)", "enzyme-lowering allele lowers (q < 0.05)"), "no evidence (q >= 0.05)")]
cat("co-effect family size:", nrow(co), "; q < 0.05:", sum(co$q < 0.05), "\n")
cw <- dcast(co, id ~ trait, value.var = c("beta_lower", "se", "q", "flag"))
g <- merge(g, cw, by.x = "vid", by.y = "id", all.x = TRUE)
for (tr in c("tg_glgc", "bmi_yengo", "homair_manning")) g[is.na(get(paste0("flag_", tr))), paste0("flag_", tr) := "variant not in file"]
g[, fasting_insulin := "not available with direction (on-disk MAGIC FI file has Bayes factors only)"]
g[, `:=`(ldl = "not downloaded", t2d = "not downloaded", cad = "not downloaded")]

# ---------------------------------------------------------------- LOEUF (gnomAD v4.1)
gn <- fread(file.path(ROOT, "data/external/gnomad_constraint/gnomad.v4.1.constraint_metrics.tsv"),
            select = c("gene", "gene_id", "transcript", "canonical", "mane_select", "lof.oe_ci.upper"))
gn[, pri := fifelse(mane_select %in% c(TRUE, "true", "TRUE"), 1L, fifelse(canonical %in% c(TRUE, "true", "TRUE"), 2L, 3L))]
gn <- gn[pri < 3 & is.finite(`lof.oe_ci.upper`)][order(pri)]
ge <- gn[grepl("^ENSG", gene_id)][!duplicated(gene_id), .(ensembl = gene_id, loeuf_e = `lof.oe_ci.upper`)]
gs <- gn[!duplicated(gene), .(symbol = gene, loeuf_s = `lof.oe_ci.upper`)]
g <- merge(g, ge, by = "ensembl", all.x = TRUE); g <- merge(g, gs, by = "symbol", all.x = TRUE)
g[, loeuf := fifelse(is.finite(loeuf_e), loeuf_e, loeuf_s)]
g[, loeuf_source := fifelse(is.finite(loeuf_e), "gnomAD v4.1 (Ensembl ID)", fifelse(is.finite(loeuf_s), "gnomAD v4.1 (symbol)", "not in gnomAD table"))]
g[, c("loeuf_e", "loeuf_s") := NULL]

# ---------------------------------------------------------------- LoF burden
hb <- h26[, .(ensembl, ALT_M1.1_beta, ALT_M1.1_p_value, AST_M1.1_beta, AST_M1.1_p_value, chen_z, chen_p,
              lineB_test, lineB_effect, lineB_q, lineB_call, lineB_status)]
g <- merge(g, hb, by = "ensembl", all.x = TRUE)
rd_burden <- function(f, beta_col) {
  b <- fread(f)
  b <- b[effect_allele == "M1.1"]
  b[, ensembl := sub("^.*\\((ENSG[0-9]+)\\).*$", "\\1", Name)]
  if (beta_col == "odds_ratio") b[, beta := log(odds_ratio)]
  b[, .(ensembl, beta, p = p_value)]
}
for (x in list(list("GGT", "data/external/GCST90083014/GCST90083014_buildGRCh38.tsv.gz", "beta"),
               list("NALD", "data/external/GCST90085487/GCST90085487_buildGRCh38.tsv.gz", "odds_ratio"))) {
  b <- rd_burden(file.path(ROOT, x[[2]]), x[[3]])
  b <- b[ensembl %in% g$ensembl & is.finite(beta) & is.finite(p)]
  b[, q := p.adjust(p, "BH")]   # family: this test over the H26 call genes with a row
  setnames(b, c("beta", "p", "q"), paste0(x[[1]], "_M1.1_", c("beta", "p", "q")))
  g <- merge(g, b, by = "ensembl", all.x = TRUE)
}
lab <- function(beta, q) fifelse(!is.finite(beta), "no burden row (too few carriers)",
                        fifelse(q < 0.05, fifelse(beta < 0, "LoF lowers trait (q < 0.05)", "LoF raises trait (q < 0.05)"),
                                fifelse(beta < 0, "LoF lowers trait (q >= 0.05)", "LoF raises trait (q >= 0.05)")))
g[, GGT_M1.1_direction := lab(GGT_M1.1_beta, GGT_M1.1_q)]
g[, NALD_M1.1_direction := lab(NALD_M1.1_beta, NALD_M1.1_q)]

# ---------------------------------------------------------------- hepatocyte share (not ambient-corrected)
hs <- fread(file.path(ROOT, "Analysis/Hypotheses/DIG1-20261024T113951Z/out_hep_share/gene_hep_share.tsv.gz"))
hs <- dcast(hs, gene ~ cohort, value.var = "median_hep_share", fun.aggregate = function(x) x[1])
setnames(hs, setdiff(names(hs), "gene"), paste0("hep_share_uncorrected_", setdiff(names(hs), "gene")))
g <- merge(g, hs, by.x = "symbol", by.y = "gene", all.x = TRUE)
g[, hep_share_ambient_corrected := "not available (no ambient-corrected per-gene table on disk)"]

# ---------------------------------------------------------------- liver protein (PXD051911)
px <- fread(file.path(ROOT, "data/PXD051911/liver_protein_quant.txt"), na.strings = c("NA", "", "NaN"))
sc <- names(px)[-(1:3)]
px[, frac := rowMeans(!is.na(as.matrix(.SD))), .SDcols = sc]
pl <- px[, .(symbol = unlist(strsplit(Genes, ";"))), by = .(ProteinAccessions, frac)][, .(liver_protein_frac_samples = max(frac)), by = symbol]
g <- merge(g, pl, by = "symbol", all.x = TRUE)
g[, liver_protein_detected := fifelse(is.finite(liver_protein_frac_samples) & liver_protein_frac_samples > 0, "yes", "not detected")]
g[is.na(liver_protein_frac_samples), liver_protein_frac_samples := 0]

# ---------------------------------------------------------------- flags and output
g[, instrument_study_non_eur := study_superpop != "EUR"]
g[, defining_trait_not_enzyme := defining_trait %in% c("NAFLD", "PDFF")]
g[, liver_fat_test_not_independent := defining_trait == "PDFF"]
setorder(g, call, combined, symbol)
fwrite(g, file.path(OUT, "dir6_gene_classes.tsv"), sep = "\t", na = "NA")

# ---------------------------------------------------------------- controls
cw2 <- w[is_ctrl]
ctl <- rbind(
  cw2[id %in% c("ctrl_PNPLA3_rs738409", "ctrl_TM6SF2_rs58542926", "ctrl_GCKR_rs1260326"),
      .(control = "concordant on liver fat expected", id, call_pdff, lfsr_pdff, pm_pdff, call_cirr, lfsr_cirr, pass = call_pdff == "concordant")],
  cw2[grepl("HSD17B13", id), .(control = "HSD17B13: ALT-raising allele concordant on cirrhosis, no liver-fat support", id, call_pdff, lfsr_pdff, pm_pdff, call_cirr, lfsr_cirr,
                               pass = call_cirr == "concordant" & call_pdff == "no outcome support")],
  cw2[grepl("GPT|GOT1", id), .(control = "enzyme-only locus: no lfsr < 0.05 on either outcome", id, call_pdff, lfsr_pdff, pm_pdff, call_cirr, lfsr_cirr,
                               pass = call_pdff == "no outcome support" & call_cirr %in% c("no outcome support", "untestable"))])
fwrite(ctl, file.path(OUT, "controls.tsv"), sep = "\t"); print(ctl)

# ---------------------------------------------------------------- calibration (descriptive)
cal <- data.table(symbol = c("PNPLA3", "HSD17B13", "MTARC1", "CIDEB", "INHBE", "GPAM"))
cal <- merge(cal, h26[, .(symbol, h26_class = primary_evidence_class, lineA_status)], by = "symbol", all.x = TRUE)
cal <- merge(cal, g[, .(symbol, call, defining_trait, combined, call_pdff, lfsr_pdff, call_cirr, lfsr_cirr, published_coding_mechanism)], by = "symbol", all.x = TRUE)
cal[is.na(combined), combined := fifelse(is.na(h26_class), "not in the H26 table (no DIR6 class)", "not among the 286 calls (no DIR6 class)")]
fwrite(cal, file.path(OUT, "calibration.tsv"), sep = "\t", na = "NA"); print(cal)

print(counts[set == "all 286 calls" & call == "all"])
writeLines(c(sprintf("bound_pdff\t%.5f", bounds["pdff"]), sprintf("bound_cirr\t%.5f", bounds["cirr"]), paste0("bound_cirr_rule\t", bound_rule_cirr),
             sprintf("pi0_pdff\t%.4f", get_pi0(fits$pdff)), sprintf("pi0_cirr\t%.4f", get_pi0(fits$cirr)),
             sprintf("independent_loci\t%d", uniqueN(g$locus_id)), sprintf("unique_instruments\t%d", uniqueN(g$vid))),
           file.path(OUT, "run_constants.tsv"))
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
cat("STEP2 DONE\n")
