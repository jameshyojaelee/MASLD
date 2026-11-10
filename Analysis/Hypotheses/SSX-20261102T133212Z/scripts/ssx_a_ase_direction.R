# SSX-A (PRESPEC_SSX_A.md; exploratory). Do GTEx liver eQTL directions hold in diseased liver? For s3h sites that are
# themselves GTEx liver significant variants for the same gene, a heterozygous donor carries the ALT allele on one copy,
# so the pooled ALT read fraction should move in the direction of the GTEx slope (no phasing needed).
# Het rule (S3h): REF + ALT >= 10 and 0.1 < ALT fraction < 0.9. Unit: included SteatoSITE donor (515).
# Site choice (no allele fractions): most het donors, then higher median het depth, then lower position.
# Testable: >= 30 het donors. Shift = qlogis(pooled ALT fraction) - median qlogis(pooled fraction) of testable null
#   sites in the same median-depth quintile x EUR MAF bin (fallback to depth quintile if < 20 null sites; S3h rule).
# Primary: share of testable genes with sign(shift) = sign(slope); one-sided p under 10,000 random sign draws.
# Usage: Rscript ssx_a_ase_direction.R --out <new dir>
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE)
opt <- function(k) { i <- match(k, a); if (is.na(i)) NA_character_ else a[i + 1] }
out <- opt("--out"); stopifnot(!is.na(out), !dir.exists(out)); dir.create(out, recursive = TRUE)
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
SS <- file.path(ROOT, "Analysis/Hypotheses/SS-20261025T135237Z")
SSX <- file.path(ROOT, "Analysis/Hypotheses/SSX-20261102T133212Z")
GT <- file.path(ROOT, "data/external/allelic_refs/gtex_v8/GTEx_Analysis_v8_eQTL")
stopifnot(system(paste("cd", SSX, "&& sha256sum -c --quiet PRESPEC_SSX.sha256")) == 0)
SEED <- 20261102L

U <- fread(file.path(SS, "out_s2_qc/units.tsv"))[toupper(as.character(include)) == "TRUE"]
U[, sample_alias := as.character(sample_alias)]
st <- fread(cmd = paste("tr -d '\\r' <", shQuote(file.path(ROOT, "data/external/SteatoSite/fastq/PRJEB58625_fibrosis_stage_20261025T133450Z.tsv"))),
            colClasses = "character")
U[, stage_num := suppressWarnings(as.numeric(st$nash_crn_kleiner_fibrosis_stage[match(sample_alias, st$sample_alias)]))]
TG <- fread(file.path(SS, "out_s0d_targets/s1b_targets_annotated.tsv"))[role %in% c("s3h", "s3h_null")]
TG <- TG[, .(chrom, pos, ref, alt, role, gene, id, af_eur, maf = pmin(af_eur, 1 - af_eur))]
TG[, variant_id := paste(chrom, pos, ref, alt, "b38", sep = "_")]
A <- fread(file.path(SS, "out_s3fgh-20261026T233706Z/counts/ase_counts.tsv.gz"), colClasses = c(sample_alias = "character"))
A <- merge(A[sample_alias %in% U$sample_alias], TG, by = c("chrom", "pos"))
A[, dp := nref + nalt][, f := nalt / dp]
A <- merge(A, U[, .(sample_alias, stage_num)], by = "sample_alias")
H <- A[dp >= 10 & f > 0.1 & f < 0.9]

# GTEx liver: significant variant-gene pairs (slope sign) and eGene leads (log2 aFC)
SG <- fread(cmd = paste("zcat", shQuote(file.path(GT, "Liver.v8.signif_variant_gene_pairs.txt.gz"))),
            select = c("variant_id", "gene_id", "slope", "pval_nominal"))
SG[, id := sub("\\..*", "", gene_id)]
EG <- fread(cmd = paste("zcat", shQuote(file.path(GT, "Liver.v8.egenes.txt.gz"))),
            select = c("gene_id", "variant_id", "log2_aFC", "qval"))
EG[, id := sub("\\..*", "", gene_id)]

# site summaries over all included donors
site_sum <- function(h) h[, .(n_het = .N, med_dp = as.numeric(median(dp)), pooled_f = sum(nalt) / sum(dp)),
                          by = .(role, id, chrom, pos, variant_id, maf)]
S_all <- site_sum(H)
S_all[, unit := fifelse(role == "s3h", id, paste0("null:", chrom, ":", pos))]

# eQTL-eligible sites: s3h site that is a GTEx liver significant variant for the same gene
EL <- merge(S_all[role == "s3h"], SG[, .(variant_id, id, slope, pval_nominal)], by = c("variant_id", "id"))
setorder(EL, id, -n_het, -med_dp, pos)
EQ <- EL[, .SD[1], by = id]
EQ[, testable := n_het >= 30]
NUL <- S_all[role == "s3h_null"][, testable := n_het >= 30]
n_pairs <- nrow(merge(unique(TG[role == "s3h", .(variant_id, id)]), SG[, .(variant_id, id)], by = c("variant_id", "id")))

# bins: depth quintile breaks from all testable sites (eQTL-chosen + null), MAF bins as S3h
bins <- function(d, brk) {
  d[, dbin := cut(med_dp, brk, include.lowest = TRUE, labels = FALSE)]
  d[, mbin := cut(maf, c(0.05, 0.15, 0.30, 0.50), include.lowest = TRUE, right = FALSE, labels = FALSE)] }
null_centre <- function(d, nl, col) {
  vapply(seq_len(nrow(d)), function(i) {
    pool <- nl[dbin == d$dbin[i] & mbin == d$mbin[i]][[col]]
    if (length(pool) < 20) pool <- nl[dbin == d$dbin[i]][[col]]
    if (!length(pool)) NA_real_ else median(pool) }, numeric(1)) }
ql <- function(p) qlogis(pmin(pmax(p, 1e-4), 1 - 1e-4))

score <- function(eq, nl) {   # eq, nl: testable sites with pooled_f, med_dp, maf
  brk <- unique(quantile(c(eq$med_dp, nl$med_dp), seq(0, 1, 0.2)))
  eq <- bins(copy(eq), brk); nl <- bins(copy(nl), brk)
  nl[, plogit := ql(pooled_f)]; eq[, plogit := ql(pooled_f)]
  eq[, null_med_logit := null_centre(eq, nl, "plogit")]
  eq[, null_med_f := null_centre(eq, nl, "pooled_f")]
  eq[, shift := plogit - null_med_logit][, agree := sign(shift) == sign(slope)]
  eq }
perm_share <- function(n, obs, B = 10000L) {
  set.seed(SEED); nullv <- rbinom(B, n, 0.5) / n
  list(p = (1 + sum(nullv >= obs)) / (1 + B), null_mean = mean(nullv), null_q95 = unname(quantile(nullv, 0.95))) }

T1 <- score(EQ[testable == TRUE], NUL[testable == TRUE])
n_t <- nrow(T1); share <- mean(T1$agree)
pp <- perm_share(n_t, share)
bt <- binom.test(sum(T1$agree), n_t, 0.5, alternative = "greater")$p.value

# beside: lead-variant aFC
T1[, is_lead := variant_id == EG$variant_id[match(id, EG$id)]]
T1[, aFC := EG$log2_aFC[match(id, EG$id)]][, exp_f := 2^aFC / (1 + 2^aFC)]
T1[, obs_f_centred := pooled_f - (null_med_f - 0.5)]
LD <- T1[is_lead == TRUE & is.finite(aFC)]
sp <- if (nrow(LD) >= 5) cor.test(LD$exp_f, LD$obs_f_centred, method = "spearman", exact = FALSE) else NULL

# beside: F0-F1 and F3-F4 donors separately (same sites; null centring recomputed within the same donor subset)
by_stage <- function(keep, label) {
  h <- H[stage_num %in% keep]; s <- site_sum(h)
  eq <- merge(T1[, .(id, chrom, pos, slope)], s[role == "s3h"], by = c("id", "chrom", "pos"))
  nl <- s[role == "s3h_null" & paste(chrom, pos) %in% paste(NUL[testable == TRUE]$chrom, NUL[testable == TRUE]$pos)]
  r <- score(eq, nl)
  rbind(data.table(subset = label, min_het = 1, n_genes = nrow(r), share = mean(r$agree, na.rm = TRUE)),
        data.table(subset = label, min_het = 10, n_genes = sum(r$n_het >= 10), share = mean(r[n_het >= 10]$agree, na.rm = TRUE)))
}
STG <- rbind(by_stage(1, "F0-F1"), by_stage(c(3, 4), "F3-F4"))

# sensitivity: RNA-based heterozygote deficit F > 0.1; median het depth above the 99th percentile
G <- A[dp >= 10, .(n = .N, n_het = sum(f > 0.1 & f < 0.9), n_alt_hom = sum(f >= 0.9)), by = .(chrom, pos)]
G[, p_alt := (2 * n_alt_hom + n_het) / (2 * n)][, F_het := 1 - (n_het / n) / (2 * p_alt * (1 - p_alt))]
T1 <- merge(T1, G[, .(chrom, pos, F_het)], by = c("chrom", "pos"), all.x = TRUE)
dp99 <- quantile(c(T1$med_dp, NUL[testable == TRUE]$med_dp), 0.99)
sens <- function(keep, label) { r <- T1[keep]; data.table(analysis = label, n_genes = nrow(r), share = mean(r$agree),
  p_perm = perm_share(nrow(r), mean(r$agree))$p) }
SEN <- rbind(sens(T1$F_het <= 0.1 | !is.finite(T1$F_het), "drop F_het > 0.1"),
             sens(T1$med_dp <= dp99, "drop median het depth > 99th percentile"),
             sens((T1$F_het <= 0.1 | !is.finite(T1$F_het)) & T1$med_dp <= dp99, "both"))

# secondary (assay check): share of het donors with binomial p < 0.05 against the bin null median fraction,
# S3h-chosen catalog sites vs null sites; label permutation within depth x MAF bins
S3 <- copy(S_all[role == "s3h"]); setorder(S3, unit, -n_het, -med_dp, pos); S3 <- S3[, .SD[1], by = unit]
SEC <- rbind(S3[n_het >= 30], NUL[testable == TRUE], fill = TRUE)[, testable := NULL]
brk2 <- unique(quantile(SEC$med_dp, seq(0, 1, 0.2))); SEC <- bins(SEC, brk2)
SEC[, null_med_f := null_centre(SEC, SEC[role == "s3h_null"], "pooled_f")]
HH <- merge(H, SEC[, .(chrom, pos, null_med_f)], by = c("chrom", "pos"))
HH[, bp := mapply(function(x, n, p) binom.test(x, n, p)$p.value, nalt, dp, null_med_f)]
SEC <- merge(SEC, HH[, .(sig_share = mean(bp < 0.05)), by = .(chrom, pos)], by = c("chrom", "pos"))
d_obs <- SEC[role == "s3h", mean(sig_share)] - SEC[role == "s3h_null", mean(sig_share)]
set.seed(SEED); SEC[, cell := paste(dbin, mbin)]
d_null <- replicate(1000, { lab <- SEC[, .(role = sample(role), sig_share), by = cell]
  lab[role == "s3h", mean(sig_share)] - lab[role == "s3h_null", mean(sig_share)] })
p_sec <- (1 + sum(d_null >= d_obs)) / 1001

# synthetic check: planted signed imbalance on the testable eQTL sites, and random signs
syn <- function(effect, random_sign) {
  set.seed(SEED + 1L)
  h <- merge(H[role == "s3h"], T1[, .(chrom, pos, slope, null_med_logit)], by = c("chrom", "pos"))
  sgn <- if (random_sign) h[, .(chrom, pos)][, s := sample(c(-1, 1), 1), by = .(chrom, pos)]$s else sign(h$slope)
  h[, nalt := rbinom(.N, dp, plogis(null_med_logit + effect * sgn))]
  s <- site_sum(h)
  eq <- merge(T1[, .(id, chrom, pos, slope)], s, by = c("id", "chrom", "pos"))
  mean(score(eq, NUL[testable == TRUE])$agree) }
SYN <- data.table(check = c("planted 0.3 logit in slope direction", "planted 0.3 logit random sign"),
                  share = c(syn(0.3, FALSE), syn(0.3, TRUE)))

fwrite(T1[order(id)], file.path(out, "primary_genes.tsv"), sep = "\t")
fwrite(EL, file.path(out, "eligible_sites.tsv"), sep = "\t")
fwrite(STG, file.path(out, "beside_by_stage.tsv"), sep = "\t")
fwrite(SEN, file.path(out, "sensitivity.tsv"), sep = "\t")
fwrite(SEC, file.path(out, "secondary_sites.tsv"), sep = "\t")
fwrite(SYN, file.path(out, "synthetic_check.tsv"), sep = "\t")
summ <- c(
  paste("included donors:", nrow(U), "; s3h-GTEx significant site-gene pairs:", n_pairs,
        "; eligible genes with >= 1 het:", uniqueN(EL$id), "; testable genes (>= 30 het donors):", n_t),
  paste("PRIMARY share with GTEx slope sign:", sprintf("%.4f", share), "(", sum(T1$agree), "of", n_t, ")",
        "; random-sign null mean", sprintf("%.4f", pp$null_mean), "95th pct", sprintf("%.4f", pp$null_q95),
        "; p perm", signif(pp$p, 3), "; exact binomial p", signif(bt, 3),
        if (n_t < 30) "; DESCRIPTIVE (fewer than 30 testable genes)" else ""),
  paste("lead-site Spearman (expected vs centred observed ALT fraction): n", nrow(LD),
        if (!is.null(sp)) paste("rho", sprintf("%.3f", sp$estimate), "p", signif(sp$p.value, 3)) else "not computed"),
  paste("secondary: mean share of hets with binomial p < 0.05, catalog", sprintf("%.4f", SEC[role == "s3h", mean(sig_share)]),
        "(n", SEC[, sum(role == "s3h")], ") vs null", sprintf("%.4f", SEC[role == "s3h_null", mean(sig_share)]),
        "(n", SEC[, sum(role == "s3h_null")], "); p perm", signif(p_sec, 3)),
  paste("null testable sites:", NUL[, sum(testable)], "; null pooled ALT fraction median:", signif(median(NUL[testable == TRUE]$pooled_f), 4)))
writeLines(summ, file.path(out, "summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
options(width = 200); cat(summ, sep = "\n"); print(STG); print(SEN); print(SYN)
