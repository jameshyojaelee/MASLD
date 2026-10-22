# H51 stage 1 (blinded). No fit in this script uses the observed GSE213621 bin labels as a design term.
# GSE213621 observed labels enter only: filterByExpr(group = bin) (smallest group size only), stratum
# and cell definitions for permutation and bootstrap draws, and label-only count tables.
source("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H51-20261015T112617Z/h51_common.R")
t0 <- Sys.time()
ts <- format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC")
OUT <- if (DEBUG) file.path(H, paste0("debug-", ts), "stage1") else file.path(H, "stage1")
if (dir.exists(OUT)) stop("Output dir exists, refusing: ", OUT)
dir.create(OUT, recursive = TRUE)
log_msg("Stage 1 start. DEBUG =", DEBUG, " OUT =", OUT, " cores =", NCORES, " B_PERM =", B_PERM, " B_BOOT =", B_BOOT)
verify_hashes(OUT, c("h51_common.R", "stage1_blinded.R", "stage2_outcome.R"))
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

X <- load_inputs()
dge_all <- X$dge_all; meta <- X$meta; ann <- X$ann; disc <- X$disc; g <- X$g
md <- character()  # diagnostics report lines
add <- function(...) md <<- c(md, paste0(...))
tbl <- function(dt) { dt <- as.data.frame(dt); c(paste0("| ", paste(names(dt), collapse = " | "), " |"),
  paste0("|", paste(rep("---", ncol(dt)), collapse = "|"), "|"),
  apply(dt, 1, function(r) paste0("| ", paste(trimws(format(r)), collapse = " | "), " |"))) }

## ---------- Diagnostic 2: counts ----------
stopifnot(nrow(disc) == 371L, nrow(g) == 258L)
nd <- disc[, .N, by = .(dataset, bin, inferred_sex)][order(dataset, bin, inferred_sex)]
n3831 <- fread(INPUTS[id == "n3831", path])
chk_n <- merge(nd, n3831, by = c("dataset", "bin", "inferred_sex"), all = TRUE, suffixes = c("", "_3831"))
stopifnot(identical(chk_n$N, chk_n$N_3831))
add("# H51 stage 1 diagnostics (blinded)\n\nRun dir: ", H, "\nOutput: ", OUT, "\nDEBUG: ", DEBUG, "\n")
add("## 2. Unit counts\n\nDiscovery disease units (unit = genetic-identity cluster), by cohort, bin, sex (identical to job 3831 n table: TRUE):\n")
add(tbl(dcast(nd, dataset + inferred_sex ~ bin, value.var = "N", fill = 0)))
add("\nGSE213621 disease units by bin and sex:\n")
add(tbl(dcast(g[, .N, by = .(inferred_sex, bin)], inferred_sex ~ bin, value.var = "N", fill = 0)))
st <- dcast(g[, .N, by = .(acc_q, inferred_sex, bin)], acc_q + inferred_sex ~ bin, value.var = "N", fill = 0)
add("\nPermutation strata (accession quartile x inferred sex), units per bin:\n"); add(tbl(st))
S4_fires <- any(as.matrix(st[, ..BINS]) == 0)
add("\nS4 (stratum lacks a bin): ", if (S4_fires) "FIRES" else "does not fire", "; minimum stratum-bin count ", min(as.matrix(st[, ..BINS])), "\n")
if (S4_fires) stop("S4 fired: halves replacement not implemented because prespec counted it as not expected; stopping for review")
# Cross-check strata against the validity file
ab <- fread(INPUTS[id == "acc_bins", path])
xc <- merge(g[, .(sample_id, acc_rank, acc_q = as.character(acc_q), bin = as.character(bin), inferred_sex)], ab, by = "sample_id", suffixes = c("", ".v"))
stopifnot(nrow(xc) == 258L, all(xc$acc_rank == xc$acc_rank.v), all(xc$acc_q == xc$acc_q.v), all(xc$bin == xc$bin.v), all(xc$inferred_sex == xc$inferred_sex.v))
# condition agrees with fibrosis_stage coding
stopifnot(all(as.character(g$bin) == c("F0F1", "F2", "F3F4")[g$fibrosis_stage]))
add("Strata and bins match `accession_order_bins.tsv` for all 258 units; `condition` agrees with stage codes 1/2/3.\n")
for (h in DISC) { s <- disc[dataset != h]; o <- disc[dataset == h]
  add("LOCO hold out ", h, ": selection ", nrow(s), " (", paste(table(s$bin), collapse = "/"), "), held out ", nrow(o), " (", paste(table(o$bin), collapse = "/"), ")") }
add("")

## ---------- Random draws (all, serially, fixed order) ----------
set.seed(SEED, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
strata <- sort(unique(g$stratum))
perm <- matrix(NA_character_, B_PERM, nrow(g), dimnames = list(NULL, g$sample_id))
for (k in seq_len(B_PERM)) {
  lab <- as.character(g$bin)
  for (s in strata) { i <- which(g$stratum == s); lab[i] <- lab[i][sample.int(length(i))] }
  perm[k, ] <- lab
}
boot_draw <- function(cells) {
  cl <- sort(unique(cells))
  unlist(lapply(cl, function(cc) { i <- which(cells == cc); i[sample.int(length(i), length(i), replace = TRUE)] }), use.names = FALSE)
}
bootG <- t(vapply(seq_len(B_BOOT), function(k) boot_draw(g$cell), integer(nrow(g))))
hsets <- lapply(setNames(DISC, DISC), function(h) disc[dataset == h])
bootH <- list()
for (h in DISC) bootH[[h]] <- t(vapply(seq_len(B_BOOT), function(k) boot_draw(hsets[[h]]$cell), integer(nrow(hsets[[h]]))))
# write
pl <- data.table(draw = seq_len(B_PERM), perm); fwrite(pl, file.path(OUT, "perm_labels.tsv.gz"), sep = "\t")
bl <- rbindlist(c(list(data.table(set = "GSE213621", draw = rep(seq_len(B_BOOT), each = nrow(g)), position = rep(seq_len(nrow(g)), B_BOOT), index = as.vector(t(bootG)))),
  lapply(DISC, function(h) data.table(set = h, draw = rep(seq_len(B_BOOT), each = nrow(hsets[[h]])), position = rep(seq_len(nrow(hsets[[h]])), B_BOOT), index = as.vector(t(bootH[[h]]))))))
fwrite(bl, file.path(OUT, "boot_indices.tsv.gz"), sep = "\t")
# check: permutation preserves stratum x bin counts; bootstrap preserves cell sizes
stopifnot(all(apply(perm, 1, function(l) all(table(g$stratum, factor(l, BINS)) == table(g$stratum, g$bin)))))
stopifnot(all(apply(bootG, 1, function(i) all(table(g$cell[i]) == table(g$cell)))))
log_msg("Random draws written: perms", B_PERM, "boot", B_BOOT)

## ---------- Wave A fits: D0, LOCO selection, held-out single-cohort ----------
gene_ids <- rownames(dge_all)
summ_fit <- function(f, coefs) {
  rbindlist(lapply(coefs, function(cf) data.table(gene = rownames(f$coefficients), contrast = cf, logFC = unname(f$coefficients[, cf]),
    SE = unname(f$stdev.unscaled[, cf] * sqrt(f$s2.post)), P = unname(f$p.value[, cf]), FDR = p.adjust(f$p.value[, cf], "BH"))))
}
jobsA <- c(list(list(type = "D0", h = NA)), lapply(DISC, function(h) list(type = "sel", h = h)), lapply(DISC, function(h) list(type = "ho", h = h)))
resA <- mclapply(jobsA, function(j) {
  if (j$type == "D0") {
    d <- subdge(dge_all, disc$sample_id); f <- fit_qw(d, des_D0(disc)); return(list(tab = summ_fit(f, CF)))
  }
  if (j$type == "sel") {
    s <- disc[dataset != j$h]; d <- subdge(dge_all, s$sample_id); f <- fit_qw(d, des_D0(s)); return(list(tab = summ_fit(f, CF)))
  }
  o <- hsets[[j$h]]; d <- subdge(dge_all, o$sample_id)
  keep <- filterByExpr(d, group = o$bin); d <- d[keep, , keep.lib.sizes = TRUE]
  f <- fit_qw(d, des_G0(o)); list(tab = summ_fit(f, CF), keep = gene_ids[keep])
}, mc.cores = min(NCORES, length(jobsA)), mc.preschedule = FALSE)
if (any(vapply(resA, inherits, TRUE, "try-error"))) { print(resA); stop("Wave A fit failed") }
D0 <- resA[[1]]$tab
SEL <- setNames(lapply(2:4, function(i) resA[[i]]$tab), DISC)
HO <- setNames(lapply(5:7, function(i) resA[[i]]$tab), DISC)
KEEP_H <- setNames(lapply(5:7, function(i) resA[[i]]$keep), DISC)
log_msg("Wave A fits done", format(Sys.time() - t0))

## ---------- S1: reproduce job 3831 ----------
s3831 <- fread(INPUTS[id == "sum3831", path])
D0s <- D0[, .(n_units = nrow(disc), n_genes = .N, n_bh05 = sum(FDR < 0.05), n_bh05_up = sum(FDR < 0.05 & logFC > 0),
              n_bh05_down = sum(FDR < 0.05 & logFC < 0), median_abs_logFC_hits = median(abs(logFC[FDR < 0.05]))), by = contrast]
s1 <- merge(D0s, s3831[, .(contrast, n_units, n_bh05, n_bh05_up, n_bh05_down, median_abs_logFC_hits)], by = "contrast", suffixes = c("", "_3831"))
s1[, ok := n_units == n_units_3831 & abs(n_bh05 - n_bh05_3831) <= 2 & abs(n_bh05_up - n_bh05_up_3831) <= 2 &
        abs(n_bh05_down - n_bh05_down_3831) <= 2 & round(median_abs_logFC_hits, 4) == round(median_abs_logFC_hits_3831, 4)]
S1_pass <- nrow(s1) == 2L && all(s1$ok)
add("## 1. S1: D0 reproduction of job 3831\n"); add(tbl(fmt(s1, 6))); add("\nS1: ", if (S1_pass) "PASS" else "FAIL", "\n")
cat(if (S1_pass) "S1 REPRODUCTION PASSED\n" else "S1 REPRODUCTION FAILED\n"); print(s1)
D0[, gene_name := ann$gene_name[match(gene, ann$gene_id_versioned)]]
D0[, chromosome := ann$chromosome[match(gene, ann$gene_id_versioned)]]
fwrite(D0, file.path(OUT, "D0_discovery_grouped_per_gene.tsv.gz"), sep = "\t")

write_stop <- function(rule, text) {
  dec <- list(stopped = TRUE, rule = rule, text = text)
  jsonlite::write_json(dec, file.path(OUT, "stage1_decision.json"), auto_unbox = TRUE, pretty = TRUE)
  writeLines(md, file.path(OUT, "stage1_diagnostics.md"))
  res <- c("# RESULT_H51", "", paste0("Verdict: stopped_by_rule (", rule, ")"), "", text, "",
           "Figures 3E and 3F remain three-cohort results; nothing about transport is claimed.")
  rf <- file.path(if (DEBUG) dirname(OUT) else H, "RESULT_H51.md")
  if (file.exists(rf)) stop("RESULT exists, refusing to overwrite: ", rf)
  writeLines(res, rf)
  writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
  quit(save = "no", status = 0)
}
if (!S1_pass) write_stop("S1", paste("D0 did not reproduce job 3831 within tolerance:", paste(capture.output(print(s1)), collapse = "\n")))

## ---------- GSE213621 gene filter, S2, S6 ----------
dgG_all <- subdge(dge_all, g$sample_id)
keepG_l <- filterByExpr(dgG_all, group = g$bin)
KEEP_G <- gene_ids[keepG_l]
dgG <- dgG_all[keepG_l, , keep.lib.sizes = TRUE]
stopifnot(identical(dgG$samples$lib.size, dgG_all$samples$lib.size))
HC <- lapply(CF, function(cf) D0[contrast == cf & FDR < 0.05, gene])
HCk <- lapply(HC, function(x) intersect(x, KEEP_G))
S2_min <- c(F3F4 = 1000L, F2 = 300L)
S2 <- data.table(contrast = names(CF), D0_hits = lengths(HC), hits_kept_G = lengths(HCk), minimum = S2_min[names(CF)])
S2[, testable := hits_kept_G >= minimum]
autos <- paste0("chr", 1:22)
nc1_genes <- HCk$F3F4[ann$chromosome[match(HCk$F3F4, ann$gene_id_versioned)] %in% autos]
add("## 3. Genes kept\n\nGSE213621 filterByExpr(group = bin), defaults: ", length(KEEP_G), " of 24,115 kept.\n")
for (h in DISC) add("Held-out ", h, " filterByExpr(group = bin): ", length(KEEP_H[[h]]), " kept.")
add("\nS2 (hits kept in GSE213621):\n"); add(tbl(S2))
add("\nNC1 gene set: F3F4 hits kept in GSE213621 on chr1-22: ", length(nc1_genes), " (removed as chrX/chrY/other: ", length(HCk$F3F4) - length(nc1_genes),
    "; of which chrX ", sum(ann$chromosome[match(HCk$F3F4, ann$gene_id_versioned)] == "chrX"), ", chrY ", sum(ann$chromosome[match(HCk$F3F4, ann$gene_id_versioned)] == "chrY"), ")\n")
PCG <- c("COL1A1", "LUM", "EPCAM", "CYP2C19", "GNMT")
pcid <- ann[gene_name %in% PCG, .(gene_name, gene = gene_id_versioned)]
stopifnot(nrow(pcid) == 5L)
pcid[, kept_G := gene %in% KEEP_G]
S6_n <- sum(pcid$kept_G); S6_pass <- S6_n >= 4L

## ---------- Null fits (stratified permutations, GSE213621) ----------
log_msg("Null fits start: ", B_PERM)
null_list <- mclapply(seq_len(B_PERM), function(k) {
  info <- data.frame(inferred_sex = g$inferred_sex, bin = factor(perm[k, ], levels = BINS))
  f <- fit_qw(dgG, des_G0(info))
  out <- list(coef = f$coefficients[, c("binF2", "binF3F4", SEXCF)])
  if (k == 1L) out$se <- cbind(binF2 = f$stdev.unscaled[, "binF2"] * sqrt(f$s2.post), binF3F4 = f$stdev.unscaled[, "binF3F4"] * sqrt(f$s2.post), df_total = f$df.total)
  out
}, mc.cores = NCORES, mc.preschedule = TRUE)
if (any(vapply(null_list, inherits, TRUE, "try-error"))) stop("Null fit failed")
NULLC <- lapply(c(binF2 = "binF2", binF3F4 = "binF3F4", sex = SEXCF), function(cf) sapply(null_list, function(z) z$coef[, cf]))
SE1 <- null_list[[1]]$se
rm(null_list); invisible(gc())
saveRDS(NULLC, file.path(OUT, "null_coefficients.rds"))
saveRDS(SE1, file.path(OUT, "null_draw1_se.rds"))
log_msg("Null fits done", format(Sys.time() - t0))

## ---------- Held-out bootstrap fits ----------
jobsB <- unlist(lapply(DISC, function(h) lapply(seq_len(B_BOOT), function(k) list(h = h, k = k))), recursive = FALSE)
bootH_list <- mclapply(jobsB, function(j) {
  o <- hsets[[j$h]]; idx <- bootH[[j$h]][j$k, ]
  ob <- o[idx]; d <- subdge(dge_all, o$sample_id, genes = KEEP_H[[j$h]])[, idx]
  colnames(d) <- make.unique(colnames(d))
  f <- fit_qw(d, des_G0(ob)); f$coefficients[, c("binF2", "binF3F4")]
}, mc.cores = NCORES, mc.preschedule = TRUE)
if (any(vapply(bootH_list, inherits, TRUE, "try-error"))) stop("Held-out bootstrap fit failed")
BOOTH <- lapply(setNames(DISC, DISC), function(h) { w <- which(vapply(jobsB, function(j) j$h == h, TRUE))
  lapply(c(binF2 = "binF2", binF3F4 = "binF3F4"), function(cf) sapply(bootH_list[w], function(m) m[, cf])) })
rm(bootH_list); invisible(gc())
saveRDS(BOOTH, file.path(OUT, "heldout_boot_coefficients.rds"))
log_msg("Held-out bootstraps done", format(Sys.time() - t0))

## ---------- Diagnostic 8: stratified null summaries, S3 ----------
x_of <- function(cf, genes) D0[contrast == cf][match(genes, gene), logFC]
nullsum <- rbindlist(lapply(names(CF), function(cn) {
  cf <- CF[[cn]]; gs <- HCk[[cn]]; x <- x_of(cf, gs); ii <- match(gs, KEEP_G)
  sl <- apply(NULLC[[cf]][ii, , drop = FALSE], 2, ols_slope, x = x)
  cc <- apply(NULLC[[cf]][ii, , drop = FALSE], 2, conc, x = x)
  data.table(contrast = cn, n_genes = length(gs), null_mean = mean(sl), null_sd = sd(sl), null_q95 = unname(quantile(sl, 0.95)),
             MDE_slope = MDE_SLOPE_K * sd(sl), conc_null_mean = mean(cc), conc_null_sd = sd(cc), n_draws = length(sl))
}))
nullsum[, S3_untestable := MDE_slope > 0.30]
nullsum[, null_mean_in_band := abs(null_mean) <= 0.05]
add("## 8. Stratified permutation null (", B_PERM, " draws; x = D0 log2FC over D0 hits kept in GSE213621)\n")
add(tbl(fmt(nullsum, 4)))
add("\nPrespecified diagnostic prediction: null mean within [-0.05, 0.05]. A mean outside is the slope that submission order and sex structure alone produce.\n")

## ---------- Diagnostic 10: LOCO reference, S5, s_c ----------
loco <- list(); loco_genes <- list()
for (cn in names(CF)) { cf <- CF[[cn]]; comp <- cn
  for (h in DISC) {
    sel <- SEL[[h]][contrast == cf]; hs <- sel[FDR < 0.05, gene]
    gs <- intersect(intersect(hs, KEEP_H[[h]]), KEEP_G)
    x <- sel[match(gs, gene), logFC]
    yh <- HO[[h]][contrast == cf][match(gs, gene), logFC]
    bh <- ols_slope(yh, x)
    bb <- apply(BOOTH[[h]][[cf]][match(gs, KEEP_H[[h]]), , drop = FALSE], 2, ols_slope, x = x)
    nb <- table(hsets[[h]]$bin)
    elig <- nb[["F0F1"]] >= 10 && nb[[comp]] >= 10
    loco_genes[[paste(cn, h)]] <- gs
    loco[[paste(cn, h)]] <- data.table(contrast = cn, heldout = h, sel_hits = length(hs), genes_used = length(gs), n_F0F1 = nb[["F0F1"]], n_comp = nb[[comp]],
      eligible = elig, b_h = bh, boot_se = sd(bb), ci_norm_lo = bh - z975 * sd(bb), ci_norm_hi = bh + z975 * sd(bb),
      ci_pct_lo = unname(quantile(bb, 0.025)), ci_pct_hi = unname(quantile(bb, 0.975)))
  } }
loco <- rbindlist(loco)
S5 <- loco[, .(n_eligible = sum(eligible), mean_bh_eligible = if (any(eligible)) mean(b_h[eligible]) else NA_real_), by = contrast]
S5[, size_determinable := n_eligible >= 2 & !is.na(mean_bh_eligible) & mean_bh_eligible >= 0.10]
S5[, s_c := fifelse(size_determinable, pmin(pmax(mean_bh_eligible, 0), 1), 1)]
add("## 10. Leave-one-cohort-out reference (b_h: held-out single-cohort log2FC on two-cohort selection log2FC)\n")
add(tbl(fmt(loco, 4)))
add("\nS5 and s_c (s_c = mean eligible b_h bounded to [0,1]; 1 if S5 fails):\n"); add(tbl(fmt(S5, 4)))

## ---------- Diagnostic 9: blinded per-gene SE and MDE ----------
d9 <- rbindlist(lapply(names(CF), function(cn) { cf <- CF[[cn]]; gs <- HCk[[cn]]; ii <- match(gs, KEEP_G)
  se <- SE1[ii, cf]; mde <- MDE_GENE_K * se; e <- S5[contrast == cn, s_c] * abs(x_of(cf, gs))
  data.table(contrast = cn, n_genes = length(gs), median_SE = median(se), median_MDE_g = median(mde), share_MDE_le_e = mean(mde <= e), s_c = S5[contrast == cn, s_c]) }))
add("## 9. Blinded per-gene SE (first null draw)\n"); add(tbl(fmt(d9, 4)))

## ---------- Diagnostic 11: positive controls feasibility, S6 ----------
pc <- merge(pcid, D0[contrast == "binF3F4", .(gene, D0_F3F4_log2FC = logFC, D0_F3F4_FDR = FDR)], by = "gene")
pc <- merge(pc, D0[contrast == "binF2", .(gene, D0_F2_log2FC = logFC)], by = "gene")
pc[, blinded_SE_F3F4 := ifelse(kept_G, SE1[match(gene, KEEP_G), "binF3F4"], NA_real_)]
pc[, blinded_MDE_nominal := MDE_SLOPE_K * blinded_SE_F3F4]
pc[, frac_units_count_ge10 := rowMeans(dgG_all$counts[gene, , drop = FALSE] >= 10)]
add("## 11. Positive-control feasibility (S6)\n"); add(tbl(fmt(pc, 4)))
add("\nS6: ", S6_n, " of 5 kept; ", if (S6_pass) "PASS" else "FAIL (primary reported as untestable)", "\n")

## ---------- Diagnostic 4: library complexity ----------
cx_units <- rbind(disc[, .(sample_id, dataset, bin)], g[, .(sample_id, dataset, bin)])
cnt <- dge_all$counts[, cx_units$sample_id]
lib <- colSums(cnt); p <- sweep(cnt, 2, lib, "/")
cx_units[, `:=`(libsize = lib, top50 = apply(p, 2, function(x) sum(sort(x, decreasing = TRUE)[1:50])),
                entropy = apply(p, 2, function(x) { x <- x[x > 0]; -sum(x * log2(x)) }), detected = colSums(cnt >= 10))]
rm(cnt, p)
cxs <- cx_units[, .(n = .N, med_lib_M = median(libsize) / 1e6, med_top50 = median(top50), med_entropy = median(entropy), med_detected = median(detected)), by = .(dataset, bin)][order(dataset, bin)]
cxp <- cx_units[, .(kruskal_p_lib = kruskal.test(libsize ~ bin)$p.value, kruskal_p_top50 = kruskal.test(top50 ~ bin)$p.value,
                    kruskal_p_entropy = kruskal.test(entropy ~ bin)$p.value, kruskal_p_detected = kruskal.test(detected ~ bin)$p.value), by = dataset]
add("## 4. Library complexity by bin (descriptive)\n"); add(tbl(fmt(cxs, 4)))
add(""); add(tbl(fmt(cxp, 3)))
fwrite(cx_units, file.path(OUT, "library_complexity_units.tsv"), sep = "\t")

## ---------- Diagnostic 5: accession order ----------
chi <- chisq.test(table(g$acc_q, g$bin))$p.value
lc <- cpm(dgG_all, log = TRUE, prior.count = 2); lc <- lc[rowMeans(lc) > 2, ]
top <- lc[order(-apply(lc, 1, var))[1:2000], ]
pcs <- prcomp(t(top), center = TRUE, scale. = FALSE)
pcr <- data.table(PC = 1:6, rho_accession = sapply(1:6, function(k) cor(pcs$x[, k], g$acc_rank, method = "spearman")),
                  p = sapply(1:6, function(k) cor.test(pcs$x[, k], g$acc_rank, method = "spearman", exact = FALSE)$p.value))
add("## 5. Accession order\n\nAccession quartile x bin:\n"); add(tbl(dcast(g[, .N, by = .(acc_q, bin)], acc_q ~ bin, value.var = "N")))
add("\nChi-square p = ", signif(chi, 3), "\n\nLabel-free expression PCs (top 2,000 variable genes, mean logCPM > 2) vs accession rank:\n")
add(tbl(fmt(pcr, 3))); rm(lc, top)

## ---------- Diagnostic 6: composition (discovery cohorts only) ----------
comp <- rbindlist(lapply(DISC, function(h) { o <- hsets[[h]]; z <- load_ilr(h, o$sample_id)
  r2 <- summary(lm(as.integer(o$bin) ~ z$ilr))$r.squared
  data.table(dataset = h, n = nrow(o), share_below_1e6_floor = z$share_below_floor, R2_ordinal_bin_on_6ILR = r2) }))
zG <- load_ilr("GSE213621", g$sample_id)
comp <- rbind(comp, data.table(dataset = "GSE213621", n = nrow(g), share_below_1e6_floor = zG$share_below_floor, R2_ordinal_bin_on_6ILR = NA_real_))
add("## 6. Composition (7-part lineage coding, floor 1e-6, 6 ILR coordinates). GSE213621 R2 is computed in stage 2 only.\n")
add(tbl(fmt(comp, 4)))

## ---------- Diagnostic 7: label reliability ----------
sra <- fread(INPUTS[id == "sra", path])
xl <- merge(g[, .(sample_id, condition)], sra[, .(Run, fibrotic_stage)], by.x = "sample_id", by.y = "Run")
lab_ok <- nrow(xl) == 258L && all(cond_to_bin[xl$condition] == xl$fibrotic_stage)
pr <- fread(INPUTS[id == "pairs", path]); setnames(pr, sub("^#", "", names(pr)))
ro <- fread(INPUTS[id == "roster", path]); ds <- setNames(ro$dataset, ro$sample_id); lab <- setNames(sra$fibrotic_stage, sra$Run)
q <- pr[ds[IID1] == "GSE213621" & ds[IID2] == "GSE213621" & KINSHIP > 0.45]
q[, `:=`(l1 = lab[IID1], l2 = lab[IID2], in1 = IID1 %in% meta$sample_id, in2 = IID2 %in% meta$sample_id)]
q[, combo := mapply(function(a, b) paste(sort(c(a, b)), collapse = "|"), l1, l2)]
add("## 7. Label reliability\n\nDisease units whose `condition` matches SRA `fibrotic_stage`: ", sum(cond_to_bin[xl$condition] == xl$fibrotic_stage), " of 258 (all: ", lab_ok, ").")
add("GSE213621 same-genotype pairs (kinship > 0.45): ", nrow(q), "; discordant bins ", sum(q$l1 != q$l2), "; discordant pairs with a member in the 774 roster ", sum(q$l1 != q$l2 & (q$in1 | q$in2)), "; pairs with both in roster ", sum(q$in1 & q$in2), ".\n")
add(tbl(q[, .N, by = .(same = l1 == l2, combo)][order(same, combo)]))

## ---------- Decision and save ----------
dec <- list(stopped = FALSE, S1_pass = S1_pass, S4_fires = S4_fires, S6_n_kept = S6_n, S6_pass = S6_pass,
            S2 = S2, S3 = nullsum[, .(contrast, MDE_slope, S3_untestable)], S5 = S5, debug = DEBUG)
jsonlite::write_json(dec, file.path(OUT, "stage1_decision.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
saveRDS(list(KEEP_G = KEEP_G, KEEP_H = KEEP_H, HC = HC, HCk = HCk, nc1_genes = nc1_genes, SEL = SEL, HO = HO,
             loco = loco, loco_genes = loco_genes, S5 = S5, S2 = S2, nullsum = nullsum, pc = pc, S6_pass = S6_pass,
             S1 = s1, comp = comp, cx_units = cx_units, d9 = d9),
        file.path(OUT, "stage1_objects.rds"))
add("\n## Stop-rule summary\n")
add("- S1 reproduction: ", if (S1_pass) "PASS" else "FAIL")
for (cn in names(CF)) add("- ", cn, ": S2 hits kept ", S2[contrast == cn, hits_kept_G], " (min ", S2[contrast == cn, minimum], ") -> ", if (S2[contrast == cn, testable]) "testable" else "UNTESTABLE",
  "; S3 MDE_slope ", signif(nullsum[contrast == cn, MDE_slope], 3), " -> ", if (nullsum[contrast == cn, S3_untestable]) "UNTESTABLE" else "testable",
  "; S5 size ", if (S5[contrast == cn, size_determinable]) "determinable" else "not determinable", " (eligible ", S5[contrast == cn, n_eligible], ", mean b_h ", signif(S5[contrast == cn, mean_bh_eligible], 3), ")")
add("- S4: ", if (S4_fires) "fires" else "does not fire")
add("- S6: ", S6_n, " of 5 positive controls kept -> ", if (S6_pass) "PASS" else "FAIL")
add("\nElapsed: ", format(Sys.time() - t0))
writeLines(md, file.path(OUT, "stage1_diagnostics.md"))
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
log_msg("STAGE 1 DONE", format(Sys.time() - t0))
