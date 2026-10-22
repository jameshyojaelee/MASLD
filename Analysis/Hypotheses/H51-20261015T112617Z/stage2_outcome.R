# H51 stage 2 (outcome). Reads stage 1 draws; draws no random numbers (except DEBUG label scrambling).
source("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H51-20261015T112617Z/h51_common.R")
t0 <- Sys.time()
S1DIR <- Sys.getenv("H51_STAGE1_DIR", file.path(H, "stage1"))
ts <- format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC")
OUT <- file.path(dirname(S1DIR), paste0("stage2-", ts))
if (dir.exists(OUT)) stop("Output dir exists, refusing: ", OUT)
dir.create(OUT)
log_msg("Stage 2 start. DEBUG =", DEBUG, " S1DIR =", S1DIR, " OUT =", OUT)
SCRIPTS <- c("h51_common.R", "stage1_blinded.R", "stage2_outcome.R")
verify_hashes(OUT, SCRIPTS)
s1s <- fread(file.path(S1DIR, "script_sha256.tsv")); s2s <- fread(file.path(OUT, "script_sha256.tsv"))
if (!identical(s1s$sha256, s2s$sha256)) stop("Scripts changed since stage 1; refusing (record in DEVIATIONS.md and rerun deliberately)")
dec1 <- jsonlite::read_json(file.path(S1DIR, "stage1_decision.json"))
if (isTRUE(dec1$stopped)) { log_msg("Stage 1 stopped by rule", dec1$rule, "; stage 2 does nothing."); quit(save = "no", status = 0) }
if (isTRUE(dec1$debug) != DEBUG) stop("DEBUG flag differs between stage 1 and stage 2")

X <- load_inputs(); dge_all <- X$dge_all; ann <- X$ann; disc <- X$disc; g <- X$g
S <- readRDS(file.path(S1DIR, "stage1_objects.rds"))
NULLC <- readRDS(file.path(S1DIR, "null_coefficients.rds"))
BOOTH <- readRDS(file.path(S1DIR, "heldout_boot_coefficients.rds"))
D0 <- fread(file.path(S1DIR, "D0_discovery_grouped_per_gene.tsv.gz"))
bl <- fread(file.path(S1DIR, "boot_indices.tsv.gz"))
bootG <- bl[set == "GSE213621"][order(draw, position)]
B <- max(bootG$draw); bootG <- matrix(bootG$index, nrow = B, byrow = TRUE)
stopifnot(ncol(bootG) == nrow(g))
KEEP_G <- S$KEEP_G; HCk <- S$HCk; HC <- S$HC

if (DEBUG) {  # never use observed labels while debugging: stratified scramble with a separate seed
  set.seed(99L); lab <- as.character(g$bin)
  for (s in sort(unique(g$stratum))) { i <- which(g$stratum == s); lab[i] <- lab[i][sample.int(length(i))] }
  g[, bin := factor(lab, levels = BINS)]; g[, cell := paste(bin, inferred_sex, sep = "_")]
  log_msg("DEBUG: GSE213621 labels replaced by a stratified scramble")
}

dgG <- subdge(dge_all, g$sample_id, genes = KEEP_G)
ilrG <- load_ilr("GSE213621", g$sample_id)$ilr
ilrD <- do.call(rbind, lapply(DISC, function(h) load_ilr(h, disc[dataset == h, sample_id])$ilr))[disc$sample_id, ]
disc[, top50z := { cnt <- dge_all$counts[, sample_id]; p <- sweep(cnt, 2, colSums(cnt), "/")
  t5 <- apply(p, 2, function(x) sum(sort(x, decreasing = TRUE)[1:50])); as.numeric(scale(t5)) }, by = dataset]

pack <- function(f, cfs) list(coef = f$coefficients[, cfs, drop = FALSE],
  se = sapply(cfs, function(cf) f$stdev.unscaled[, cf] * sqrt(f$s2.post)), p = f$p.value[, cfs, drop = FALSE], df = f$df.total)
CFS3 <- c("binF2", "binF3F4", SEXCF)
jobs <- c(list("G0", "G0_acc", "G0_comp", "G0_rle", "G0_female", "D0_top50", "D0_no162694", "D0_comp"), as.list(paste0("boot", seq_len(B))))
res <- mclapply(jobs, function(j) {
  if (j == "G0") return(pack(fit_qw(dgG, des_G0(g)), CFS3))
  if (j == "G0_acc") return(pack(fit_qw(dgG, model.matrix(~ factor(inferred_sex) + acc_q + bin, g)), CFS3))
  if (j == "G0_comp") { info <- data.frame(g[, .(inferred_sex, bin)], ilrG); return(pack(fit_qw(dgG, model.matrix(~ factor(inferred_sex) + ilr1 + ilr2 + ilr3 + ilr4 + ilr5 + ilr6 + bin, info)), CFS3)) }
  if (j == "G0_rle") { d <- calcNormFactors(dgG, method = "RLE"); return(pack(fit_qw(d, des_G0(g)), CFS3)) }
  if (j == "G0_female") { w <- which(g$inferred_sex == "F"); return(pack(fit_qw(dgG[, w], model.matrix(~ bin, g[w])), c("binF2", "binF3F4"))) }
  if (j == "D0_top50") return(pack(fit_qw(subdge(dge_all, disc$sample_id), model.matrix(~ factor(dataset) + factor(inferred_sex) + top50z + bin, disc)), c("binF2", "binF3F4")))
  if (j == "D0_no162694") { s <- disc[dataset != "GSE162694"]; return(pack(fit_qw(subdge(dge_all, s$sample_id), des_D0(s)), c("binF2", "binF3F4"))) }
  if (j == "D0_comp") { info <- data.frame(disc[, .(dataset, inferred_sex, bin)], ilrD)
    return(pack(fit_qw(subdge(dge_all, disc$sample_id), model.matrix(~ factor(dataset) + factor(inferred_sex) + ilr1 + ilr2 + ilr3 + ilr4 + ilr5 + ilr6 + bin, info)), c("binF2", "binF3F4"))) }
  k <- as.integer(sub("boot", "", j)); idx <- bootG[k, ]; d <- dgG[, idx]; colnames(d) <- make.unique(colnames(d))
  fit_qw(d, des_G0(g[idx]))$coefficients[, CFS3]
}, mc.cores = NCORES, mc.preschedule = FALSE)
if (any(vapply(res, inherits, TRUE, "try-error"))) { print(res[vapply(res, inherits, TRUE, "try-error")]); stop("Stage 2 fit failed") }
names(res) <- unlist(jobs)
FIT <- res[1:8]
BOOTC <- lapply(setNames(CFS3, c("binF2", "binF3F4", "sex")), function(cf) sapply(res[-(1:8)], function(m) m[, cf]))
rm(res); invisible(gc())
saveRDS(FIT, file.path(OUT, "observed_fits.rds")); saveRDS(BOOTC, file.path(OUT, "G0_boot_coefficients.rds"))
log_msg("Fits done", format(Sys.time() - t0))

## ---------- helpers ----------
xD0 <- function(cf, gs) D0[contrast == cf][match(gs, gene), logFC]
yv <- function(fit, cf, gs) fit$coef[match(gs, rownames(fit$coef)), cf]
iG <- function(gs) match(gs, KEEP_G)
slope_block <- function(y, x, gs, cfnull, cfboot) {
  b <- ols_slope(y, x)
  nl <- if (!is.null(cfnull)) apply(NULLC[[cfnull]][iG(gs), , drop = FALSE], 2, ols_slope, x = x) else NA
  bb <- apply(BOOTC[[cfboot]][iG(gs), , drop = FALSE], 2, ols_slope, x = x)
  data.table(n_genes = length(gs), slope = b, null_mean = mean(nl), null_sd = sd(nl),
             p_perm_onesided = if (!is.null(cfnull)) (1 + sum(nl >= b)) / (length(nl) + 1) else NA_real_,
             boot_se = sd(bb), ci_norm_lo = b - z975 * sd(bb), ci_norm_hi = b + z975 * sd(bb),
             ci_pct_lo = unname(quantile(bb, 0.025)), ci_pct_hi = unname(quantile(bb, 0.975)))
}
tbl <- function(dt) { dt <- as.data.frame(dt); c(paste0("| ", paste(names(dt), collapse = " | "), " |"),
  paste0("|", paste(rep("---", ncol(dt)), collapse = "|"), "|"),
  apply(dt, 1, function(r) paste0("| ", paste(trimws(format(r)), collapse = " | "), " |"))) }

## ---------- Primary ----------
S2 <- as.data.table(S$S2); NS <- as.data.table(S$nullsum); S5 <- as.data.table(S$S5)
prim <- rbindlist(lapply(names(CF), function(cn) { cf <- CF[[cn]]; gs <- HCk[[cn]]; x <- xD0(cf, gs); y <- yv(FIT$G0, cf, gs)
  nl_c <- apply(NULLC[[cf]][iG(gs), , drop = FALSE], 2, conc, x = x)
  cbind(data.table(contrast = cn), slope_block(y, x, gs, cf, cf), conc_obs = conc(y, x), conc_null_mean = mean(nl_c), conc_null_sd = sd(nl_c)) }))
prim[, testable := S2[match(contrast, S2$contrast), testable] & !NS[match(contrast, NS$contrast), S3_untestable]]
prim[, p_holm := p.adjust(ifelse(testable, p_perm_onesided, NA_real_), "holm", n = 2)]
prim[, holm_pass := testable & !is.na(p_holm) & p_holm < 0.05]

## ---------- Positive controls (S7) ----------
pc <- as.data.table(S$pc)
pc[, `:=`(G0_F3F4_log2FC = ifelse(kept_G, FIT$G0$coef[match(gene, KEEP_G), "binF3F4"], NA_real_),
          G0_F3F4_p = ifelse(kept_G, FIT$G0$p[match(gene, KEEP_G), "binF3F4"], NA_real_),
          G0_F2_log2FC = ifelse(kept_G, FIT$G0$coef[match(gene, KEEP_G), "binF2"], NA_real_))]
pc[, sign_ok := sign(G0_F3F4_log2FC) == sign(D0_F3F4_log2FC)]
pc[, sign_and_p := kept_G & sign_ok & G0_F3F4_p < 0.05]
pc[, ratio_to_D0 := abs(G0_F3F4_log2FC) / abs(D0_F3F4_log2FC)]
S7_n <- sum(pc$sign_and_p, na.rm = TRUE); S7_pass <- S7_n >= 3L
pred_all_sign <- all(pc$sign_ok[pc$kept_G]); pred_ge4_p <- sum(pc$kept_G & pc$G0_F3F4_p < 0.05, na.rm = TRUE) >= 4
pred_ratio <- all(pc$ratio_to_D0[pc$kept_G] >= 0.25 & pc$ratio_to_D0[pc$kept_G] <= 1.5)
f2pred <- pc[gene_name %in% c("COL1A1", "LUM", "CYP2C19") & kept_G]
f2_ok <- all(sign(f2pred$G0_F2_log2FC) == c(COL1A1 = 1, LUM = 1, CYP2C19 = -1)[f2pred$gene_name] & abs(f2pred$G0_F2_log2FC) < abs(f2pred$G0_F3F4_log2FC))

## ---------- NC1 ----------
nc <- S$nc1_genes; xnc <- xD0("binF3F4", nc)
NC1 <- slope_block(yv(FIT$G0, SEXCF, nc), xnc, nc, NULL, "sex")
NC1[, fails := (ci_norm_lo > 0 | ci_norm_hi < 0) & abs(slope) > 0.10]

## ---------- Size reference R_c ----------
loco <- as.data.table(S$loco)
lf <- rbindlist(lapply(seq_len(nrow(loco)), function(r) { L <- loco[r]; cf <- CF[[L$contrast]]; gs <- S$loco_genes[[paste(L$contrast, L$heldout)]]
  sel <- S$SEL[[L$heldout]][contrast == cf]; x <- sel[match(gs, gene), logFC]
  bG <- ols_slope(yv(FIT$G0, cf, gs), x); bGb <- apply(BOOTC[[cf]][iG(gs), , drop = FALSE], 2, ols_slope, x = x)
  bhb <- apply(BOOTH[[L$heldout]][[cf]][match(gs, S$KEEP_H[[L$heldout]]), , drop = FALSE], 2, ols_slope, x = x)
  data.table(contrast = L$contrast, heldout = L$heldout, eligible = L$eligible, genes = length(gs), b_h = L$b_h, b_G = bG,
    b_G_lo = bG - z975 * sd(bGb), b_G_hi = bG + z975 * sd(bGb), ratio = bG / L$b_h, ratio_boot_se = sd(bGb / bhb),
    diff = bG - L$b_h, diff_lo = (bG - L$b_h) - z975 * sd(bGb - bhb), diff_hi = (bG - L$b_h) + z975 * sd(bGb - bhb),
    bGb = list(bGb), bhb = list(bhb)) }))
Rc <- rbindlist(lapply(names(CF), function(cn) { e <- lf[contrast == cn & eligible == TRUE]
  if (nrow(e) == 0) return(data.table(contrast = cn, n_eligible = 0L))
  num <- Reduce(`+`, e$bGb) / nrow(e); den <- Reduce(`+`, e$bhb) / nrow(e); Rb <- num / den
  R <- mean(e$b_G) / mean(e$b_h)
  data.table(contrast = cn, n_eligible = nrow(e), mean_b_G = mean(e$b_G), mean_b_h = mean(e$b_h), R = R, boot_se = sd(Rb),
             ci_norm_lo = R - z975 * sd(Rb), ci_norm_hi = R + z975 * sd(Rb), ci_pct_lo = unname(quantile(Rb, 0.025)),
             ci_pct_hi = unname(quantile(Rb, 0.975)), share_den_below_0.05 = mean(den < 0.05)) }), fill = TRUE)
Rc[, size_determinable := S5[match(contrast, S5$contrast), size_determinable]]
Rc[, size_reading := fifelse(!size_determinable, "not determinable (S5)",
  fifelse(ci_norm_lo >= 0.67 & ci_norm_hi <= 1.5, "similar size", fifelse(ci_norm_hi < 1, "smaller in GSE213621",
  fifelse(ci_norm_lo > 1, "larger in GSE213621", "size not determined"))))]

## ---------- Descriptives ----------
desc <- list()
for (cn in names(CF)) { cf <- CF[[cn]]; gs <- HCk[[cn]]; x <- xD0(cf, gs)
  for (v in c("G0_acc", "G0_rle", "G0_female")) desc[[length(desc) + 1]] <- data.table(contrast = cn, analysis = v, n_genes = length(gs), slope = ols_slope(yv(FIT[[v]], cf, gs), x))
  up <- gs[x > 0]; dn <- gs[x < 0]
  desc[[length(desc) + 1]] <- cbind(data.table(contrast = cn, analysis = "up hits"), slope_block(yv(FIT$G0, cf, up), xD0(cf, up), up, cf, cf))
  desc[[length(desc) + 1]] <- cbind(data.table(contrast = cn, analysis = "down hits"), slope_block(yv(FIT$G0, cf, dn), xD0(cf, dn), dn, cf, cf))
  allg <- KEEP_G; xa <- xD0(cf, allg); ya <- yv(FIT$G0, cf, allg)
  desc[[length(desc) + 1]] <- data.table(contrast = cn, analysis = "all kept genes (selection-free)", n_genes = length(allg), slope = ols_slope(ya, xa), spearman = cor(xa, ya, method = "spearman"))
  for (v in c("D0_top50", "D0_no162694")) { hs <- intersect(rownames(FIT[[v]]$coef)[p.adjust(FIT[[v]]$p[, cf], "BH") < 0.05], KEEP_G)
    xv <- FIT[[v]]$coef[match(hs, rownames(FIT[[v]]$coef)), cf]
    desc[[length(desc) + 1]] <- cbind(data.table(contrast = cn, analysis = paste0("G0 on hits of ", v)), slope_block(yv(FIT$G0, cf, hs), xv, hs, cf, cf)) }
  hs <- intersect(rownames(FIT$D0_comp$coef)[p.adjust(FIT$D0_comp$p[, cf], "BH") < 0.05], KEEP_G)
  desc[[length(desc) + 1]] <- data.table(contrast = cn, analysis = "composition-adjusted (G0_comp on D0_comp hits)", n_genes = length(hs),
                                         slope = ols_slope(yv(FIT$G0_comp, cf, hs), FIT$D0_comp$coef[match(hs, rownames(FIT$D0_comp$coef)), cf]))
}
se3 <- fread(INPUTS[id == "stage_ext", path], select = c("gene_id_versioned", "axis", "contrast", "logFC", "FDR"))[axis == "fibrosis"]
for (pcn in c("F3_vs_F0", "F4_vs_F0")) { s <- se3[contrast == pcn & FDR < 0.05 & gene_id_versioned %in% KEEP_G]
  desc[[length(desc) + 1]] <- cbind(data.table(contrast = "F3F4", analysis = paste0("G0 F3F4 on printed Fig 3E ", pcn, " hits")),
                                    slope_block(yv(FIT$G0, "binF3F4", s$gene_id_versioned), s$logFC, s$gene_id_versioned, "binF3F4", "binF3F4")) }
desc <- rbindlist(desc, fill = TRUE)
acc_half <- merge(prim[, .(contrast, slope)], desc[analysis == "G0_acc", .(contrast, acc_slope = slope)], by = "contrast")
acc_half[, acc_below_half := acc_slope < 0.5 * slope]
comp_half <- merge(prim[, .(contrast, slope)], desc[grepl("^composition", analysis), .(contrast, comp_slope = slope)], by = "contrast")
comp_half[, comp_below_half := comp_slope < 0.5 * slope]
r2G <- summary(lm(as.integer(g$bin) ~ ilrG))$r.squared

## ---------- Per-gene table ----------
pg_state <- function(fit, cf, gs, x, sc) {
  y <- yv(fit, cf, gs); se <- fit$se[match(gs, rownames(fit$coef)), cf]; p <- fit$p[match(gs, rownames(fit$coef)), cf]
  q <- p.adjust(p, "BH"); mde <- MDE_GENE_K * se; e <- sc * abs(x)
  st <- fifelse(q < 0.05 & sign(y) == sign(x), "replicated", fifelse(q < 0.05, "opposite sign", fifelse(mde <= e, "not replicated at adequate power", "underpowered")))
  list(y = y, se = se, p = p, q = q, mde = mde, e = e, state = st, df = fit$df[match(gs, rownames(fit$coef))])
}
pgs <- rbindlist(lapply(names(CF), function(cn) { cf <- CF[[cn]]; all <- HC[[cn]]; gs <- HCk[[cn]]; sc <- S5[contrast == cn, s_c]
  x <- xD0(cf, gs); a <- pg_state(FIT$G0, cf, gs, x, sc); b <- pg_state(FIT$G0_acc, cf, gs, x, sc)
  kept <- data.table(gene = gs, gene_name = ann$gene_name[match(gs, ann$gene_id_versioned)], contrast = cn, D0_log2FC = x, G0_log2FC = a$y,
    G0_lo = a$y - qt(0.975, a$df) * a$se, G0_hi = a$y + qt(0.975, a$df) * a$se, G0_SE = a$se, G0_p = a$p, G0_BH_q = a$q, MDE_g = a$mde, e_g = a$e, s_c = sc,
    state = a$state, acc_log2FC = b$y, acc_BH_q = b$q, acc_state = b$state)
  miss <- setdiff(all, gs)
  rbind(kept, data.table(gene = miss, gene_name = ann$gene_name[match(miss, ann$gene_id_versioned)], contrast = cn, D0_log2FC = xD0(cf, miss),
    state = "not measurable in GSE213621", acc_state = "not measurable in GSE213621"), fill = TRUE) }))
pgs[, state_changed_acc := state != acc_state]
cat_ok <- setNames(vapply(names(CF), function(cn) isTRUE(prim[contrast == cn, holm_pass]) && S7_pass && isTRUE(S$S6_pass) && !isTRUE(acc_half[contrast == cn, acc_below_half]), TRUE), names(CF))
pgs[, catalog := ifelse(cat_ok[contrast], "candidate (descriptive; needs user approval)", "not proposed for the Catalog")]
if (!any(S5$size_determinable)) setnames(pgs, "e_g", "e_g_unshrunk_discovery_effect")
fwrite(pgs, file.path(OUT, "per_gene_states.tsv.gz"), sep = "\t")
stc <- pgs[, .N, by = .(contrast, state)][order(contrast, state)]
stchg <- pgs[state != "not measurable in GSE213621", .(n = .N, changed_by_accession_adjustment = sum(state_changed_acc)), by = contrast]

## ---------- Verdict (computed by the prespec rule) ----------
S6 <- S$S6_pass
verdict_c <- prim[, .(contrast, testable, p_perm_onesided, p_holm, holm_pass)]
if (!S6) { verdict <- "untestable"; rule <- sprintf("S6: %d of 5 positive controls kept by filterByExpr (minimum 4)", sum(pc$kept_G))
} else if (!any(prim$testable)) { verdict <- "untestable"; rule <- "S2/S3: no contrast testable"
} else if (!S7_pass) { verdict <- "untestable"; rule <- sprintf("S7: %d of %d kept positive controls had the discovery sign with nominal p < 0.05 (minimum 3)", S7_n, sum(pc$kept_G))
} else if (any(prim$holm_pass)) { verdict <- "pass"; rule <- paste("Holm-adjusted p < 0.05 for:", paste(prim[holm_pass == TRUE, contrast], collapse = ", "))
} else { verdict <- "null"; rule <- "No testable contrast passed Holm at family-wise 0.05; S1 to S7 passed" }

## ---------- RESULT ----------
r <- c("# RESULT_H51", "", if (DEBUG) "**DEBUG RUN: GSE213621 labels scrambled; numbers are not results.**\n" else NULL,
  paste0("Verdict: **", verdict, "**. ", rule, "."), "",
  "Unit: genetic-identity clusters (371 discovery disease units in GSE130970, GSE135251, GSE162694; 258 GSE213621 disease units). Genes are the axis of the slope, not replicates. Cross-sectional bins in unrelated biopsies.",
  "Models: limma voomWithQualityWeights on stored counts and RLE factors (774-unit). D0 `~ dataset + inferred_sex + bin`; G0 `~ inferred_sex + bin`.", "",
  "## Primary (family: 2 one-sided stratified-permutation tests, Holm, FWER 0.05)", "",
  "b = OLS slope of G0 log2FC on D0 log2FC over D0 BH hits kept in GSE213621. Null: 1,000 bin permutations within 8 accession-quartile x sex strata. Interval: 200-draw bin x sex stratified bootstrap (normal primary; percentile beside).", "",
  tbl(fmt(prim)), "", "## Stop rules", "",
  paste0("- S1 D0 reproduction of job 3831: PASS"),
  paste0("- S2 hits kept: ", paste(S2$contrast, S2$hits_kept_G, collapse = "; ")),
  paste0("- S3 MDE_slope: ", paste(NS$contrast, signif(NS$MDE_slope, 3), collapse = "; "), " (untestable if > 0.30)"),
  paste0("- S5 size reading eligible: ", paste(S5$contrast, S5$n_eligible, "cohorts, mean b_h", signif(S5$mean_bh_eligible, 3), "->", ifelse(S5$size_determinable, "determinable", "not determinable"), collapse = "; ")),
  paste0("- S6 positive controls kept: ", sum(pc$kept_G), " of 5"),
  paste0("- S7 positive controls with discovery sign and nominal p < 0.05 (G0 F3F4): ", S7_n, " -> ", if (S7_pass) "PASS" else "FAIL"), "",
  "## Positive controls (G0 F3F4 vs F0F1)", "", tbl(fmt(pc[, .(gene_name, kept_G, D0_F3F4_log2FC, G0_F3F4_log2FC, G0_F3F4_p, ratio_to_D0, G0_F2_log2FC, sign_and_p)])), "",
  paste0("Predictions: all kept have discovery sign: ", pred_all_sign, "; at least 4 with nominal p < 0.05: ", pred_ge4_p, "; each |log2FC| within 0.25-1.5 x D0: ", pred_ratio, "; F2 (COL1A1, LUM up, CYP2C19 down, smaller than F3F4): ", f2_ok), "",
  "## Negative control NC1 (G0 male-vs-female coefficient on D0 F3F4 log2FC, autosomal F3F4 hits)", "", tbl(fmt(NC1[, .(n_genes, slope, boot_se, ci_norm_lo, ci_norm_hi, fails)])),
  if (isTRUE(NC1$fails)) "NC1 failed: a GSE213621 coefficient unrelated to bin also aligns with discovery stage effects. Either sex and bin are not separated in GSE213621, or sex and fibrosis share real biology (for example sex-biased hepatocyte metabolic genes). See the female-only fit." else "NC1 passed (prediction held).", "",
  paste0("Null diagnostic (prediction: null mean within [-0.05, 0.05]): ", paste(NS$contrast, signif(NS$null_mean, 3), ifelse(NS$null_mean_in_band, "inside", "OUTSIDE"), collapse = "; ")), "",
  "## Size relative to within-discovery transport (R_c = mean b_G(-h) / mean b_h over eligible held-out cohorts)", "",
  tbl(fmt(Rc)), "", "Per fold:", "", tbl(fmt(lf[, .(contrast, heldout, eligible, genes, b_h, b_G, b_G_lo, b_G_hi, ratio, ratio_boot_se, diff, diff_lo, diff_hi)])), "",
  "## Descriptive (outside both families)", "", tbl(fmt(desc)), "",
  tbl(fmt(acc_half)), "", tbl(fmt(comp_half)), "",
  paste0("Composition: R2 of ordinal bin on the 6 ILR coordinates in GSE213621 = ", signif(r2G, 3), ". Composition-adjusted results depend on the deconvolution input (corrected_774_20261009T141331Z); no gene is named from adjusted fits."), "",
  "## Per-gene states (family: kept H_c genes per contrast, BH; reported as a family only for a contrast passing Holm)", "", tbl(stc), "", tbl(stchg), "",
  "Per-gene table: per_gene_states.tsv.gz (D0 and G0 log2FC, 95% interval, BH q, MDE_g = 3.65 x SE_g, e_g = s_c x |x_g|, state, accession-adjusted state).", "",
  "## Reading", "")
txt <- character()
for (cn in names(CF)) { P <- prim[contrast == cn]; R <- Rc[contrast == cn]
  txt <- c(txt, sprintf("- %s vs F0F1: b = %.3f (95%% bootstrap interval %.3f to %.3f), one-sided stratified-permutation p = %.4g, Holm p = %s; null mean %.3f; sign concordance %.3f (null %.3f). Size: %s (R = %.3f, %.3f to %.3f).",
    cn, P$slope, P$ci_norm_lo, P$ci_norm_hi, P$p_perm_onesided, ifelse(is.na(P$p_holm), "NA (untestable)", sprintf("%.4g", P$p_holm)), P$null_mean, P$conc_obs, P$conc_null_mean,
    R$size_reading, R$R, R$ci_norm_lo, R$ci_norm_hi)) }
if (verdict == "pass") txt <- c(txt, paste0("Grouped fibrosis-stage effects estimated in three cohorts (371 disease units) point the same way in the 258 GSE213621 disease units, beyond the slope that submission order and sex structure produce, for: ", paste(prim[holm_pass == TRUE, contrast], collapse = ", "), "."))
if (verdict == "null") txt <- c(txt, "Stage effects estimated in three cohorts do not transport detectably to the GSE213621 bins. Slopes above the upper bootstrap bound are not supported. The stage 1 MDE is the design's sensitivity (80% power), not an exclusion bound. Figures 3E and 3F remain three-cohort results. No Catalog column is proposed.")
if (verdict == "untestable") txt <- c(txt, paste0("Untestable by rule (", rule, "). Figures 3E and 3F remain three-cohort results; nothing about transport is claimed. The primary numbers above are printed only."))
for (cn in acc_half[acc_below_half == TRUE, contrast]) txt <- c(txt, paste0("- ", cn, ": accession-adjusted slope is below half the primary; the replication cannot be separated from submission order, and no Catalog column is proposed."))
for (cn in comp_half[comp_below_half == TRUE, contrast]) txt <- c(txt, paste0("- ", cn, ": composition-adjusted slope is below half the primary; the replication is carried mainly by a shared composition shift (depends on the deconvolution input)."))
txt <- c(txt, "", "Not independent support for the 11,130 stage-supported genes (RESULTS:1278), which already include GSE213621. Figure 3E legend: 3E contrasts use 387 units from three cohorts, not the 774-unit substrate (change needed regardless of H51).",
  "", paste0("Stage 1 dir: ", S1DIR, "; stage 2 dir: ", OUT, "; elapsed stage 2: ", format(Sys.time() - t0)))
r <- c(r, txt)
writeLines(r, file.path(OUT, "RESULT_H51.md"))
jsonlite::write_json(list(verdict = verdict, rule = rule, primary = prim, Rc = Rc[, !c("size_determinable")], NC1 = NC1, S7_n = S7_n, S7_pass = S7_pass, debug = DEBUG),
                     file.path(OUT, "result.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
rf <- file.path(dirname(S1DIR), "RESULT_H51.md")
if (file.exists(rf)) log_msg("RESULT_H51.md already exists at", rf, "; not overwritten") else file.copy(file.path(OUT, "RESULT_H51.md"), rf)
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
log_msg("STAGE 2 DONE. Verdict:", verdict, format(Sys.time() - t0))
