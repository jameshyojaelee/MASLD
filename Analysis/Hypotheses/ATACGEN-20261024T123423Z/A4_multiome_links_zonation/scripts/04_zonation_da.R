#!/usr/bin/env Rscript
# A4 step 4 (Q3, Q4 inputs, Q4b). See PRESPEC.md.
#   S2: paired RNA logFC (PC vs PP bin) of the 11 held-out zonation markers.
#   Q3: zone DA, PC vs PP bins of 39 donors, limma-voom ~ donor + zone, BH over tested peaks.
#   Q4 sets: ALL_HEP = DA-tested peaks; PC / PP = top 10% by +t / -t; secondary q < 0.05
#     sets; 20 seeded nulls matched on GC quintile x mean logCPM decile x width tertile.
#   Q4b: 38 donors with steatosis; ~ donor + zone + zone:steatotic (>= 5%), voom weights
#     from ~ donor + zone reused; S = #peaks with interaction p < 0.01; donor-label
#     permutations; per-bin models ~ steatotic + sex + age + well; continuous sensitivity.
# Seeds: 20261024 + 11000 + d (null draw d), + 13000 + i (permutation i).
# Usage: 04_zonation_da.R <A4_DIR> <N_PERM> <N_CORES>
suppressPackageStartupMessages({
  library(data.table); library(edgeR); library(limma); library(jsonlite)
  library(GenomicRanges); library(parallel)
})
args <- commandArgs(TRUE)
A4 <- args[1]; N_PERM <- as.integer(args[2]); NCORE <- as.integer(args[3])
SEED <- 20261024L; N_NULL <- 20L; TOP_FRAC <- 0.10; MIN_SEC <- 1000L; Q <- 0.05; P_GLOBAL <- 0.01
ATAC_MIN_COUNT <- 5; ATAC_MIN_FRAC <- 0.5; STEAT_CUT <- 5; PROM <- 1000L
SCORE_PC <- c("CYP2E1", "GLUL", "CYP1A2", "OAT"); SCORE_PP <- c("HAL", "SDS", "CPS1", "PCK1", "ASS1")
HELD_PC <- c("LGR5", "TBX3", "AHR", "SLC2A2", "ACLY", "FASN"); HELD_PP <- c("MTHFS", "GATM", "ID1", "TBX15", "SOX4")
ROOT <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
ENDPOINTS <- file.path(ROOT, "Analysis/MASLD_Model_Benchmark/executions/gse296875-phenotype-endpoints-20260825/endpoint_lock/donor_endpoints.tsv")
HG19 <- file.path(ROOT, "Analysis/Hypotheses/H22-20261009T015913Z/peaks/peaks_hg19.tsv")
W <- file.path(A4, "work"); Z <- file.path(W, "zone"); RES <- file.path(A4, "results")
LD <- file.path(W, "ldsc"); dir.create(LD, showWarnings = FALSE, recursive = TRUE)
tic <- function(...) message(format(Sys.time(), "%H:%M:%S"), " ", ...)

read_f32 <- function(path, dims) matrix(readBin(path, "numeric", n = prod(dims), size = 4, endian = "little"), nrow = dims[1])
units <- fread(file.path(Z, "units.tsv"), colClasses = list(character = "donor_id"))
dims <- fromJSON(file.path(Z, "dims.json"))
atac <- read_f32(file.path(Z, "atac_peak_by_unit.f32"), dims$atac)
rna <- read_f32(file.path(Z, "rna_gene_by_unit.f32"), dims$rna)
genes <- fread(file.path(W, "genes_tss.tsv")); pk <- fread(file.path(W, "peaks_annot.tsv"))
ep <- fread(ENDPOINTS, colClasses = list(character = "donor_id"))
stopifnot(nrow(genes) == nrow(rna), nrow(pk) == nrow(atac))
summ <- list(nuclei_binned = sum(units$nuclei), min_nuclei_unit = min(units$nuclei))

# ---- S2: held-out marker validation ----
u2 <- units[zone %in% c("PP", "PC")]
u2[, zone := factor(zone, levels = c("PP", "PC"))]; u2[, donor := factor(donor_id)]
setorder(u2, donor, zone)
c2 <- u2$unit_index + 1L
Xd <- model.matrix(~ donor + zone, data = u2)
kg <- filterByExpr(DGEList(rna[, c2]), design = Xd)
nf <- calcNormFactors(DGEList(rna[kg, c2]))$samples$norm.factors
lc <- cpm(DGEList(rna[, c2], lib.size = colSums(rna[kg, c2]), norm.factors = nf), log = TRUE, prior.count = 2)
mk <- rbind(data.table(symbol = HELD_PC, role = "held_out", expected = "PC"),
            data.table(symbol = HELD_PP, role = "held_out", expected = "PP"),
            data.table(symbol = SCORE_PC, role = "scoring", expected = "PC"),
            data.table(symbol = SCORE_PP, role = "scoring", expected = "PP"))
ipc <- which(u2$zone == "PC"); ipp <- which(u2$zone == "PP")
stopifnot(identical(as.character(u2$donor[ipc]), as.character(u2$donor[ipp])))
mk_res <- rbindlist(lapply(seq_len(nrow(mk)), function(i) {
  gi <- which(genes$symbol == mk$symbol[i])
  if (length(gi) != 1L) return(data.table(mk[i], n_symbol_matches = length(gi)))
  d <- lc[gi, ipc] - lc[gi, ipp]
  tt <- tryCatch(t.test(d), error = function(e) list(statistic = NA_real_, p.value = NA_real_))
  data.table(mk[i], n_symbol_matches = 1L, filter_by_expr = kg[gi], mean_logcpm = mean(lc[gi, c(ipc, ipp)]),
             paired_logFC = mean(d), paired_t = unname(tt$statistic), paired_p = tt$p.value,
             donors_PC_higher = sum(d > 0), donors = length(d))
}), fill = TRUE)
mk_res[, expected_sign := fifelse(expected == "PC", paired_logFC > 0, paired_logFC < 0)]
fwrite(mk_res, file.path(RES, "zone_marker_validation.tsv"), sep = "\t")
n_ok <- mk_res[role == "held_out", sum(expected_sign, na.rm = TRUE)]
summ$s2_heldout_expected_sign <- n_ok; summ$s2_heldout_total <- nrow(mk_res[role == "held_out"])
summ$s2_status <- if (n_ok >= 8) "validated" else if (n_ok >= 7) "borderline" else "not_validated"
tic("S2 held-out markers with expected sign: ", n_ok, "/11 (", summ$s2_status, ")")

# ---- Q3: zone DA ----
kp <- which(pk$primary & rowSums(atac[, c2] >= ATAC_MIN_COUNT) >= ATAC_MIN_FRAC * length(c2))
yz <- calcNormFactors(DGEList(atac[kp, c2]))
v <- voom(yz, Xd)
fit <- eBayes(lmFit(v, Xd))
tt <- topTable(fit, coef = "zonePC", number = Inf, sort.by = "none")
da <- data.table(peak_index = pk$peak_index[kp], peak_id = pk$peak_id[kp], chr = pk$chr[kp],
                 start = pk$start[kp], end = pk$end[kp], width = pk$width[kp], gc = pk$gc[kp],
                 logFC = tt$logFC, AveExpr = tt$AveExpr, t = tt$t, p = tt$P.Value, q = tt$adj.P.Val)
fwrite(da, file.path(RES, "zone_da.tsv.gz"), sep = "\t")
summ$da_tested_peaks <- nrow(da); summ$da_q05 <- da[q < Q, .N]
summ$da_q05_PC_higher <- da[q < Q & logFC > 0, .N]; summ$da_q05_PP_higher <- da[q < Q & logFC < 0, .N]
summ$da_df_residual <- fit$df.residual[1]
dep <- u2[, .(frag_per_nucleus = atac_in_peak_fragments / nuclei, umi_per_nucleus = rna_umis / nuclei), by = .(donor, zone)]
dw <- dcast(dep, donor ~ zone, value.var = c("frag_per_nucleus", "umi_per_nucleus"))
summ$median_ratio_frag_per_nucleus_PC_over_PP <- median(dw$frag_per_nucleus_PC / dw$frag_per_nucleus_PP)
summ$median_ratio_umi_per_nucleus_PC_over_PP <- median(dw$umi_per_nucleus_PC / dw$umi_per_nucleus_PP)
tic("DA tested ", nrow(da), ", q<0.05 ", summ$da_q05, " (PC higher ", summ$da_q05_PC_higher, ")")

# promoter peaks of markers (descriptive; scoring markers circular)
mg <- genes[symbol %in% mk$symbol & has_tss == TRUE]
po <- findOverlaps(GRanges(mg$chr, IRanges(mg$tss - PROM, mg$tss + PROM)), GRanges(da$chr, IRanges(da$start, da$end)))
prom <- data.table(symbol = mg$symbol[queryHits(po)], da[subjectHits(po), .(peak_id, logFC, t, q)])
prom <- merge(mk[, .(symbol, role, expected)], prom, by = "symbol", all.x = TRUE)
fwrite(prom, file.path(RES, "zone_marker_promoter_peaks.tsv"), sep = "\t")

# ---- Q4 peak sets (hg19) and matched nulls ----
h19 <- fread(HG19)
da[, lifted := peak_id %in% h19$peak_id]
n_top <- floor(TOP_FRAC * nrow(da))
o_pc <- order(-da$t)[seq_len(n_top)]; o_pp <- order(da$t)[seq_len(n_top)]
sets <- list(ALL_HEP = seq_len(nrow(da)), PC = o_pc, PP = o_pp,
             PCq = which(da$q < Q & da$logFC > 0), PPq = which(da$q < Q & da$logFC < 0))
qd <- function(x, k) cut(x, unique(quantile(x, 0:k / k)), include.lowest = TRUE, labels = FALSE)
da[, stratum := .GRP, by = .(qd(gc, 5), qd(AveExpr, 10), qd(width, 3))]
by_s <- split(seq_len(nrow(da)), da$stratum)
n_pc_s <- table(factor(da$stratum[o_pc], levels = names(by_s)))
n_pp_s <- table(factor(da$stratum[o_pp], levels = names(by_s)))
for (d in seq_len(N_NULL)) {
  set.seed(SEED + 11000L + d)
  a <- integer(0); b <- integer(0)
  for (s in names(by_s)) {
    pool <- by_s[[s]]
    if (n_pc_s[[s]]) { x <- pool[sample.int(length(pool), n_pc_s[[s]])]; a <- c(a, x); pool <- setdiff(pool, x) }
    if (n_pp_s[[s]]) b <- c(b, pool[sample.int(length(pool), n_pp_s[[s]])])
  }
  sets[[paste0("NULL_PC_", d)]] <- sort(a); sets[[paste0("NULL_PP_", d)]] <- sort(b)
}
ann <- data.table(annotation = names(sets),
                  type = sub("_[0-9]+$", "", names(sets)),
                  draw = suppressWarnings(as.integer(sub("^NULL_P[CP]_", "", names(sets)))),
                  n_peaks = lengths(sets), n_lifted = vapply(sets, function(i) sum(da$lifted[i]), 1L))
ann[, testable := TRUE]
ann[annotation %in% c("PCq", "PPq"), testable := all(n_lifted >= MIN_SEC)]
fwrite(ann, file.path(LD, "annotations.tsv"), sep = "\t")
memb <- h19[match(da$peak_id[da$lifted], peak_id), .(peak_id, chr, start, end)]
for (nm in names(sets)) set(memb, j = nm, value = as.integer(da$lifted & seq_len(nrow(da)) %in% sets[[nm]])[da$lifted])
fwrite(memb, file.path(LD, "sets.tsv.gz"), sep = "\t")
cov_chk <- rbindlist(lapply(c("PC", "PP", "NULL_PC_1", "NULL_PP_1"), function(nm)
  da[sets[[nm]], .(set = nm, mean_gc = mean(gc), mean_AveExpr = mean(AveExpr), median_width = median(width))]))
fwrite(cov_chk, file.path(RES, "zone_ldsc_set_covariate_balance.tsv"), sep = "\t")
summ$ldsc_top_n <- n_top; summ$ldsc_secondary_testable <- ann[annotation == "PCq", testable]
tic("LDSC sets written: top ", n_top, " peaks each; secondary testable ", summ$ldsc_secondary_testable)

# ---- Q4b: steatosis x zone ----
ep[, steat := suppressWarnings(as.numeric(steatosis_numeric))]
ok_d <- ep[!is.na(steat), donor_id]
u4 <- u2[donor_id %in% ok_d]; u4[, donor := droplevels(donor)]
u4 <- merge(u4, ep[, .(donor_id, steat, age = as.numeric(age_in_yr), sex = reported_sex)], by = "donor_id", sort = FALSE)
setorder(u4, donor, zone)
c4 <- u4$unit_index + 1L
X4 <- model.matrix(~ donor + zone, data = u4)
y4 <- calcNormFactors(DGEList(atac[kp, c4]))
v4 <- voom(y4, X4)
isPC <- as.numeric(u4$zone == "PC")
dlev <- levels(u4$donor)
sbin <- tapply(u4$steat, u4$donor, `[`, 1) >= STEAT_CUT
scon <- log2(1 + tapply(u4$steat, u4$donor, `[`, 1))
inter_fit <- function(lab) {
  X <- cbind(X4, zone_steat = isPC * unname(lab[as.character(u4$donor)]))
  f <- eBayes(lmFit(v4$E, X, weights = v4$weights))
  topTable(f, coef = "zone_steat", number = Inf, sort.by = "none")
}
ti <- inter_fit(sbin); S_obs <- sum(ti$P.Value < P_GLOBAL)
tc <- inter_fit(scon)
fwrite(data.table(peak_id = da$peak_id, logFC_inter = ti$logFC, t_inter = ti$t, p_inter = ti$P.Value, q_inter = ti$adj.P.Val,
                  logFC_inter_cont = tc$logFC, p_inter_cont = tc$P.Value, q_inter_cont = tc$adj.P.Val),
       file.path(RES, "zone_steatosis_interaction.tsv.gz"), sep = "\t")
tic("interaction: S_obs ", S_obs, ", q<0.05 ", sum(ti$adj.P.Val < Q), "; running ", N_PERM, " permutations")
S_perm <- unlist(mclapply(seq_len(N_PERM), function(i) {
  set.seed(SEED + 13000L + i)
  lab <- setNames(sample(unname(sbin)), names(sbin))
  sum(inter_fit(lab)$P.Value < P_GLOBAL)
}, mc.cores = NCORE, mc.preschedule = TRUE))
fwrite(data.table(perm = seq_len(N_PERM), S = S_perm), file.path(RES, "zone_steatosis_permutation_S.tsv"), sep = "\t")
summ$q4b_donors <- length(dlev); summ$q4b_steatotic_donors <- sum(sbin)
summ$q4b_S_obs <- S_obs; summ$q4b_S_perm_median <- median(S_perm)
summ$q4b_S_perm_p <- (1 + sum(S_perm >= S_obs)) / (1 + N_PERM)
summ$q4b_inter_q05 <- sum(ti$adj.P.Val < Q); summ$q4b_inter_q05_cont <- sum(tc$adj.P.Val < Q)
summ$q4b_S_obs_cont <- sum(tc$P.Value < P_GLOBAL)

per_bin <- function(zn, cont) {
  ub <- u4[zone == zn]
  ub[, x := if (cont) log2(1 + steat) else as.numeric(steat >= STEAT_CUT)]
  Xb <- model.matrix(~ x + sex + age + factor(well_id), data = ub)
  vb <- voom(calcNormFactors(DGEList(atac[kp, ub$unit_index + 1L])), Xb)
  topTable(eBayes(lmFit(vb, Xb)), coef = "x", number = Inf, sort.by = "none")
}
pb <- list(PC = per_bin("PC", FALSE), PP = per_bin("PP", FALSE), PC_cont = per_bin("PC", TRUE), PP_cont = per_bin("PP", TRUE))
fwrite(data.table(peak_id = da$peak_id, logFC_PC = pb$PC$logFC, p_PC = pb$PC$P.Value, q_PC = pb$PC$adj.P.Val,
                  logFC_PP = pb$PP$logFC, p_PP = pb$PP$P.Value, q_PP = pb$PP$adj.P.Val,
                  q_PC_cont = pb$PC_cont$adj.P.Val, q_PP_cont = pb$PP_cont$adj.P.Val),
       file.path(RES, "zone_steatosis_per_bin.tsv.gz"), sep = "\t")
for (nm in names(pb)) {
  summ[[paste0("q4b_bin_", nm, "_q05")]] <- sum(pb[[nm]]$adj.P.Val < Q)
  summ[[paste0("q4b_bin_", nm, "_p01")]] <- sum(pb[[nm]]$P.Value < P_GLOBAL)
}
summ$q4b_bin_logFC_cor_PC_PP <- cor(pb$PC$logFC, pb$PP$logFC)
write_json(summ, file.path(RES, "zone_summary.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
writeLines(capture.output(sessionInfo()), file.path(A4, "sessionInfo_04.txt"))
tic("done")
