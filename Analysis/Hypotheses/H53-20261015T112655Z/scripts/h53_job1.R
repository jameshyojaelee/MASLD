#!/usr/bin/env Rscript
# H53 job 1 (prespec "Compute", job 1): hash checks, selection (hashed before any
# caQTL row is read), GWAS lookups, GWAS-side SuSiE refits with the copied r2 runner,
# S1, then PC1/PC2/PC2b/PC3 on the three control peaks only, and S2.
# Usage: Rscript h53_job1.R [--debug-select-only]
source("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H53-20261015T112655Z/scripts/h53_common.R")
args <- commandArgs(trailingOnly = TRUE)
SELECT_ONLY <- "--debug-select-only" %in% args
NCORE <- as.integer(Sys.getenv("H53_WORKERS", "6"))
DIAG <- file.path(RUN, "diagnostics"); RES <- file.path(RUN, "results")
REFIT <- file.path(RUN, "refit")
dir.create(file.path(REFIT, "coloc"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(REFIT, "fits"), showWarnings = FALSE)
dir.create(file.path(REFIT, "logs"), showWarnings = FALSE)

stamp("H53 job 1 start; SLURM_JOB_ID =", Sys.getenv("SLURM_JOB_ID"), "; workers =", NCORE)
cat("coloc", as.character(packageVersion("coloc")), "| susieR", as.character(packageVersion("susieR")), "\n")
cat("Config: WIN", WIN, "MAF_MIN", MAF_MIN, "W", W_PRIOR, "p1", P1, "p2", P2, "p12", P12, "(alt", P12_ALT, ") MIN_SNPS", MIN_SNPS,
    "overlap.min", OVERLAP_MIN, "entry PP.H4 >", PP_ENTRY, "and LD r2 >=", LD_ENTRY, "seed", SEED, "\n")

# ---- 1. hashes ---------------------------------------------------------------
verify_inputs(cores = 16)
cat("CHECKPOINT input-hash checks PASSED\n")

# ---- 2. selection -------------------------------------------------------------
a <- load_arms(); L <- a$L; P <- a$P
pl_all <- unique(P$hit1)
po <- pl_all[!near_any(pl_all, L$hit1)]
Pp <- P[hit1 %in% po]
stopifnot(nrow(L) == 708, uniqueN(L$ensembl) == 404, nrow(unique(L[, .(gwas, hit1)])) == 581,
          length(pl_all) == 256, length(po) == 153)

# study-leads with their refit row (gene with largest r2 PP.H4 for that (study, hit1))
pick_row <- function(D) D[order(gwas, hit1, -PP.H4, gene)][, .SD[1], by = .(gwas, hit1)]
SLl <- pick_row(L)[, arm := "liver"]
SLp <- pick_row(Pp)[, arm := "platelet"]
SL <- rbind(SLl, SLp)
SL[, max_r2_pp4 := PP.H4]
setnames(SL, c("gene", "ensembl", "idx1", "idx2", "PP.H4", "hit2"), c("refit_gene", "refit_ensembl", "idx1", "idx2", "r2_pp4", "r2_hit2"))
SL[, pos19 := as.integer(sub(".*:", "", hit1))]

# GWAS lookups: p at hit1 (smallest over rows at that position) and n variants p < 5e-8 within 500 kb
studies <- unique(SL$gwas)
reg <- fread(file.path(FM, "config/gwas_registry.tsv"))
look <- rbindlist(mclapply(studies, function(st) {
  f <- file.path(FM, reg[study_name == st, sumstats_path])
  g <- fread(f, select = c("chromosome", "position", "pval"), showProgress = FALSE)
  g[, chromosome := as.character(chromosome)]
  if (!is.numeric(g$pval)) {   # some files store p below double range as text (e.g. 1e-400); parse as numeric, underflow -> 0 -> capped at 300
    pv_chr <- g$pval; g[, pval := suppressWarnings(as.numeric(pval))]
    cat("study", st, ": pval column read as text;", sum(is.na(g$pval) & !is.na(pv_chr) & pv_chr != "NA"), "values not parseable\n")
  }
  sub <- SL[gwas == st]
  rbindlist(lapply(seq_len(nrow(sub)), function(i) {
    ch <- sub$chr[i]; ps <- sub$pos19[i]
    gc <- g[chromosome == ch]
    at <- gc[position == ps, pval]
    data.table(gwas = st, hit1 = sub$hit1[i],
               gwas_p_hit1 = if (length(at)) min(at) else NA_real_,
               n_gwas_rows_at_hit1 = length(at),
               n_gws_window = gc[abs(position - ps) <= WIN & pval < 5e-8, .N])
  }))
}, mc.cores = 4, mc.preschedule = FALSE))
SL <- merge(SL, look, by = c("gwas", "hit1"), all.x = TRUE)
SL[, mlog10_lead_p := pmin(300, -log10(pmax(gwas_p_hit1, 1e-320)))]
SL[gwas_p_hit1 == 0, mlog10_lead_p := 300]

# loci
dl <- cluster_leads(SLl$hit1)
dl[, has_plt := near_any(hit1, pl_all)]
lp <- dl[, .(chr = chr[1], has_plt = any(has_plt), n_hits = .N), by = locus]
lp[, `:=`(confirmatory = !chr %in% PILOT_CHR)]
dp <- cluster_leads(po)
SL[arm == "liver", locus := paste0("L", dl$locus[match(hit1, dl$hit1)])]
SL[arm == "platelet", locus := paste0("P", dp$locus[match(hit1, dp$hit1)])]
SL[, confirmatory := !chr %in% PILOT_CHR]
SL[arm == "liver", locus_has_plt := lp$has_plt[match(sub("^L", "", locus), lp$locus)]]
SL[arm == "platelet", locus_has_plt := FALSE]
stopifnot(uniqueN(SL[arm == "liver"]$locus) == 279, uniqueN(SL[arm == "liver" & confirmatory]$locus) == 238,
          uniqueN(SL[arm == "liver" & confirmatory & !locus_has_plt]$locus) == 185,
          uniqueN(SL[arm == "platelet"]$locus) == 127, uniqueN(SL[arm == "platelet" & confirmatory]$locus) == 101)

# one lead per locus: smallest GWAS p at hit1; ties larger r2 eQTL PP.H4; then study name
choose_lead <- function(D) D[order(locus, is.na(gwas_p_hit1), gwas_p_hit1, -max_r2_pp4, gwas)][, .SD[1], by = locus]
LOC <- choose_lead(SL)
LOC[, t1_set := fifelse(arm == "liver" & confirmatory & !locus_has_plt, "T1_liver_only",
                 fifelse(arm == "platelet" & confirmatory, "T1_platelet_only",
                 fifelse(arm == "liver" & confirmatory & locus_has_plt, "shared_liver_descriptive",
                 "post_pilot_chr1_chr19")))]
SL[, selected_lead := paste(gwas, hit1) %in% paste(LOC$gwas, LOC$hit1)]

# same-biobank sensitivity: MVP_ALT_EUR loci (cluster MVP_ALT_EUR leads only)
ma <- SL[arm == "liver" & gwas == "MVP_ALT_EUR"]
dm <- cluster_leads(ma$hit1); dm[, has_plt := near_any(hit1, pl_all)]
ma[, mvp_locus := paste0("M", dm$locus[match(hit1, dm$hit1)])]
mloc <- dm[, .(chr = chr[1], has_plt = any(has_plt)), by = locus]
ma[, mvp_locus_has_plt := mloc$has_plt[match(sub("^M", "", mvp_locus), mloc$locus)]]
MVP <- ma[order(mvp_locus, is.na(gwas_p_hit1), gwas_p_hit1, -max_r2_pp4)][, .SD[1], by = mvp_locus]
MVP[, mvp_set := fifelse(confirmatory & !mvp_locus_has_plt, "MVP_ALT_liver_only_confirmatory", "other")]
stopifnot(uniqueN(ma$mvp_locus) == 87, nrow(MVP[confirmatory == TRUE]) == 75, nrow(MVP[mvp_set == "MVP_ALT_liver_only_confirmatory"]) == 57)

fwrite(SL, file.path(RES, "selection_study_leads.tsv"), sep = "\t")
fwrite(LOC, file.path(RES, "selection_loci.tsv"), sep = "\t")
fwrite(MVP[, .(mvp_locus, gwas, hit1, chr, gwas_p_hit1, mvp_set)], file.path(RES, "selection_mvp_alt_loci.tsv"), sep = "\t")
sel_sha <- sha256(file.path(RES, c("selection_study_leads.tsv", "selection_loci.tsv", "selection_mvp_alt_loci.tsv")), 1)
writeLines(paste(sel_sha, c("results/selection_study_leads.tsv", "results/selection_loci.tsv", "results/selection_mvp_alt_loci.tsv"), sep = "  "),
           file.path(DIAG, "selection_sha256.txt"))
stamp("selection written and hashed BEFORE any caQTL row is read:"); cat(readLines(file.path(DIAG, "selection_sha256.txt")), sep = "\n")

cat("\n== Diagnostics: selected leads per arm / study / chromosome ==\n")
print(LOC[, .N, by = .(t1_set)])
print(dcast(LOC[, .N, by = .(t1_set, gwas)], gwas ~ t1_set, value.var = "N", fill = 0))
print(dcast(LOC[, .N, by = .(arm, chr)], chr ~ arm, value.var = "N", fill = 0))
cat("loci on chr1/chr19 (post-pilot):", LOC[chr %in% PILOT_CHR, .N], "\n")
cat("shared confirmatory liver loci (descriptive) by study / chr:\n")
print(LOC[t1_set == "shared_liver_descriptive", .N, by = gwas]); print(LOC[t1_set == "shared_liver_descriptive", .N, by = chr][order(as.integer(chr))])
cat("selected leads with no GWAS row at hit1:", LOC[is.na(gwas_p_hit1), .N], "| study-leads:", SL[is.na(gwas_p_hit1), .N], "\n")
SL[, hit1_pos38 := lift19to38(chr[1], pos19), by = chr]
cat("hit1 liftover failures by arm:\n"); print(SL[, .(n = .N, lift_fail = sum(is.na(hit1_pos38))), by = arm])
fwrite(SL, file.path(RES, "selection_study_leads_with_lift.tsv"), sep = "\t")
if (SELECT_ONLY) { stamp("--debug-select-only: stopping"); quit(save = "no") }

# ---- 3. refits with the copied r2 runner -------------------------------------
gf <- unique(SL[, .(gwas_name = gwas, chr = as.integer(chr), ensembl = refit_ensembl)])
fwrite(gf, file.path(REFIT, "gene_filter.tsv"), sep = "\t")
tasks <- unique(gf[, .(gwas_name, chr)])[order(chr, gwas_name)]
stamp("refit tasks (study, chr):", nrow(tasks), "| (study, gene) refits:", nrow(gf))
RUNNER <- file.path(RUN, "scripts/06_susie_coloc_h53.R")
RSCRIPT <- "/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"
run_task <- function(i) {
  st <- tasks$gwas_name[i]; ch <- tasks$chr[i]
  o1 <- file.path(REFIT, "coloc", st, paste0("susie_coloc_chr", ch, ".csv"))
  o2 <- file.path(REFIT, "coloc", st, paste0("susie_coloc_pairs_chr", ch, ".tsv"))
  if (file.exists(o1) && file.exists(o2)) return(0L)
  lg <- file.path(REFIT, "logs", paste0(st, "_chr", ch, ".log"))
  env <- c(paste0("GWAS_NAME=", st), paste0("CHR_FILTER=", ch),
           paste0("COLOC_OUT_ROOT=", file.path(REFIT, "coloc")),
           paste0("COLOC_GENE_FILTER=", file.path(REFIT, "gene_filter.tsv")),
           paste0("H53_FIT_DIR=", file.path(REFIT, "fits")),
           "LD_PANEL=polyfun", paste0("EQTL_SUSIE_DIR=", EQTL_SUSIE_DIR),
           "OPENBLAS_NUM_THREADS=1", "OMP_NUM_THREADS=1", "MKL_NUM_THREADS=1")
  system2("env", c(env, RSCRIPT, RUNNER), stdout = lg, stderr = lg)
}
t0 <- Sys.time()
rc <- unlist(mclapply(seq_len(nrow(tasks)), run_task, mc.cores = NCORE, mc.preschedule = FALSE))
done <- function(i) file.exists(file.path(REFIT, "coloc", tasks$gwas_name[i], paste0("susie_coloc_pairs_chr", tasks$chr[i], ".tsv")))
fail <- which(!vapply(seq_len(nrow(tasks)), done, logical(1)))
stamp("refit pass 1 done in", format(Sys.time() - t0), "| failed tasks:", length(fail))
if (length(fail)) {
  print(tasks[fail])
  for (i in fail) run_task(i)   # sequential retry (memory)
  fail <- which(!vapply(seq_len(nrow(tasks)), done, logical(1)))
  stamp("after sequential retry, failed tasks:", length(fail)); if (length(fail)) print(tasks[fail])
}

# large inputs actually touched: AF-augmented sumstats used by the runner, eQTL fits, LD .bim
afs <- fread(file.path(FM, "config/gwas_af_sources.tsv"))
used_af <- unique(sub("/MASLD_library_design/", "/masld_atlas/", afs[study_name %in% studies, af_sumstats_path]))
eq_rds <- file.path(EQTL_SUSIE_DIR, paste0("chr", gf$chr), paste0(gf$ensembl, "_susie.rds"))
lg_all <- unlist(lapply(list.files(file.path(REFIT, "logs"), full.names = TRUE), readLines))
big <- data.table(role = c(rep("af_sumstats", length(used_af)), rep("eqtl_susie_fit", length(eq_rds)), "gwas_af_sources_config"),
                  path = c(used_af, eq_rds, file.path(FM, "config/gwas_af_sources.tsv")))
big[, exists := file.exists(path)]
big[exists == TRUE, sha256 := sha256(path, 16)]
fwrite(big, file.path(DIAG, "large_inputs_sha256.tsv"), sep = "\t")
cat("Runner log lines naming AF-augmented sumstats:", sum(grepl("Using AF-augmented", lg_all)), "\n")

# ---- 4. S1 reproduction -----------------------------------------------------
pf <- list.files(file.path(REFIT, "coloc"), pattern = "^susie_coloc_pairs_chr[0-9]+\\.tsv$", recursive = TRUE, full.names = TRUE)
RP <- rbindlist(lapply(pf, fread), fill = TRUE)
RP[, chr := as.character(chr)]
setnames(RP, c("hit1", "PP.H4"), c("refit_hit1", "refit_pp4"))
SL <- merge(SL, RP[, .(gwas, refit_ensembl = ensembl, idx1, idx2, refit_hit1, refit_pp4, refit_eligible = eligible)],
            by = c("gwas", "refit_ensembl", "idx1", "idx2"), all.x = TRUE)
SL[, fit_file := file.path(REFIT, "fits", paste0(gwas, "__", refit_ensembl, ".rds"))]
SL[, fit_exists := file.exists(fit_file)]
SL[, reproduced := fit_exists & !is.na(refit_pp4) & refit_hit1 == hit1 & abs(refit_pp4 - r2_pp4) <= 0.01]
SL[is.na(reproduced), reproduced := FALSE]
# fit diagnostics: n credible sets, top-alpha vs hit1
fd <- rbindlist(mclapply(unique(SL[fit_exists == TRUE, fit_file]), function(f) {
  x <- readRDS(f)
  data.table(fit_file = f, n_cs = length(x$cs_index), n_snps_fit = length(x$snp_ids),
             top_alpha = list(x$snp_ids[apply(x$alpha, 1, which.max)]))
}, mc.cores = 8))
SL <- merge(SL, fd[, .(fit_file, n_cs, n_snps_fit)], by = "fit_file", all.x = TRUE)
SL[fit_exists == TRUE, top_alpha_idx1 := mapply(function(f, i) fd$top_alpha[[match(f, fd$fit_file)]][i], fit_file, idx1)]
SL[, hit1_is_top_alpha := hit1 == top_alpha_idx1]
SL[, gwas_lead_snp := fifelse(!is.na(top_alpha_idx1), top_alpha_idx1, hit1)]
fwrite(SL, file.path(RES, "study_leads_job1.tsv"), sep = "\t")

LOC <- merge(LOC[, .(locus, gwas, hit1, t1_set)], SL, by = c("locus", "gwas", "hit1"), all.x = TRUE)
fwrite(LOC, file.path(RES, "loci_job1.tsv"), sep = "\t")
cat("\n== S1 reproduction ==\n")
print(SL[, .(n = .N, reproduced = sum(reproduced), rate = round(mean(reproduced), 4)), by = .(arm)])
print(SL[, .(n = .N, reproduced = sum(reproduced)), by = .(arm, gwas)][order(arm, gwas)])
cat("2022_36402844_PDFF_EUR leads (sumstats modified after r2):\n")
print(SL[gwas == "2022_36402844_PDFF_EUR", .(hit1, refit_gene, r2_pp4, refit_hit1, refit_pp4, reproduced)])
cat("Non-reproduced study-leads:\n")
print(SL[reproduced == FALSE, .(arm, gwas, hit1, refit_gene, idx1, idx2, r2_pp4, refit_hit1, refit_pp4, fit_exists)])
s1 <- LOC[t1_set %in% c("T1_liver_only", "T1_platelet_only"), .(n = .N, reproduced = sum(reproduced), rate = mean(reproduced)), by = t1_set]
print(s1)
S1_PASS <- nrow(s1) == 2 && all(s1$rate >= 0.90)
cat("S1", if (S1_PASS) "PASS" else "FAIL (stop rule)", "\n")
cat("credible sets per refit (median, IQR) by arm:\n"); print(SL[fit_exists == TRUE, as.list(quantile(n_cs, c(.25, .5, .75))), by = arm])
cat("r2 hit1 differs from top-alpha SNP of idx1, by arm:\n"); print(SL[fit_exists == TRUE, .(n = .N, differs = sum(!hit1_is_top_alpha)), by = arm])
# region SNP liftover failures by arm
lf <- rbindlist(mclapply(unique(SL[fit_exists == TRUE, .(fit_file, arm)])$fit_file, function(f) {
  gk <- gwas_keys(readRDS(f)); data.table(fit_file = f, n = nrow(gk), lift_fail = sum(is.na(gk$pos38)), non_snv = sum(!gk$snv))
}, mc.cores = 8))
lf <- merge(lf, unique(SL[, .(fit_file, arm)]), by = "fit_file")
cat("region SNP liftover failures / non-SNV by arm (summed over refits):\n"); print(lf[, .(snps = sum(n), lift_fail = sum(lift_fail), non_snv = sum(non_snv)), by = arm])
stamp("CHECKPOINT S1 reproduction computed:", if (S1_PASS) "PASS" else "FAIL")

# ---- 5. positive controls (three control peaks only) ------------------------
CTRL <- data.table(peak = c("peak21014", "peak156441", "peak273749"), chr = c("1", "19", "5"))
read_peak_rows <- function(ch, pk) {
  f <- list.files(CQDIR, pattern = paste0("^caQTL_hg38_chr", ch, "_.*txt.gz$"), full.names = TRUE)
  fread(cmd = paste0("zcat ", shQuote(f), " | awk -F'\\t' 'NR==1 || $2==\"", pk, "\"'"), sep = "\t")
}
QC <- lapply(seq_len(nrow(CTRL)), function(i) read_peak_rows(CTRL$chr[i], CTRL$peak[i])); names(QC) <- CTRL$peak
cat("caQTL rows read (control peaks only):", vapply(QC, nrow, 1L), "\n")
allq <- rbindlist(QC)
sc <- allq[MAF >= MAF_MIN, sqrt(varbeta * 2 * N_CAQTL * MAF * (1 - MAF))]
cat("sdY scale check median sqrt(varbeta*2*138*MAF*(1-MAF)) over control-peak rows (MAF>=0.01):", round(median(sc), 4),
    "IQR", round(quantile(sc, .25), 4), "-", round(quantile(sc, .75), 4), "\n")

# PC1 Currin route, pilot-coded: UKBB_GGT vs caQTL, SNVs, MAF >= 0.01, +-500 kb of the caQTL lead, sdY = 1
ggt <- fread(file.path(SSDIR, "UKBB_GGT_reformatted_hg19.tsv"), showProgress = FALSE)
ggt[, chromosome := as.character(chromosome)]
pc1 <- rbindlist(lapply(seq_len(nrow(CTRL)), function(i) {
  ch <- CTRL$chr[i]; pk <- CTRL$peak[i]
  g <- ggt[chromosome == ch & nchar(allele1) == 1 & nchar(allele2) == 1 & se > 0]
  g[, pos38 := lift19to38(ch, position)]; g <- g[!is.na(pos38)]
  g[, akey := paste(pos38, pmin(allele1, allele2), pmax(allele1, allele2), sep = ":")]; g <- g[!duplicated(akey)]
  q <- copy(QC[[pk]]); sp <- tstrsplit(q$variant, ":", fixed = TRUE)
  q[, `:=`(pos38 = as.integer(sp[[2]]), ref = sp[[3]], alt = sp[[4]])]
  q <- q[MAF >= 0.01 & nchar(ref) == 1 & nchar(alt) == 1]
  q[, akey := paste(pos38, pmin(ref, alt), pmax(ref, alt), sep = ":")]
  lead <- q[which.min(pvalue)]
  rq <- q[abs(pos38 - lead$pos38) <= WIN]; rg <- g[abs(pos38 - lead$pos38) <= WIN]
  m <- merge(rg[, .(akey, gb = beta, gv = se^2)], rq[, .(akey, qb = beta, qv = varbeta, MAF)], by = "akey")
  m <- m[is.finite(gv) & is.finite(qv) & qv > 0 & gv > 0]
  d1 <- list(beta = m$gb, varbeta = m$gv, snp = m$akey, type = "quant", sdY = 1, N = 343850)
  d2 <- list(beta = m$qb, varbeta = m$qv, snp = m$akey, type = "quant", sdY = 1, N = N_CAQTL, MAF = m$MAF)
  r5 <- suppressMessages(coloc.abf(d1, d2))                  # default p12 = 1e-5 (pilot)
  r6 <- suppressMessages(coloc.abf(d1, d2, p12 = 5e-6))
  # PC3: hand-coded lbf vs coloc lABF.df2
  hand <- caqtl_lbf(m$qb, m$qv)
  rr <- as.data.table(r5$results)
  pc3 <- max(abs(hand[match(rr$snp, m$akey)] - rr$lABF.df2))
  data.table(peak = pk, chr = ch, caqtl_lead = lead$variant, caqtl_lead_p = lead$pvalue, n_snps = nrow(m),
             PP.H3_p12_1e5 = r5$summary[["PP.H3.abf"]], PP.H4_p12_1e5 = r5$summary[["PP.H4.abf"]],
             PP.H4_p12_5e6 = r6$summary[["PP.H4.abf"]], pc3_max_abs_diff = pc3)
}))
rm(ggt); gc(FALSE)
cat("\n== PC1 (Currin route) and PC3 ==\n"); print(pc1)

# PC2 / PC2b primary route at the named peak
prim_ctrl <- function(st, h, pk) {
  r <- SL[gwas == st & hit1 == h]
  if (!nrow(r) || !isTRUE(r$fit_exists)) return(data.table(gwas = st, hit1 = h, peak = pk, note = "no refit"))
  fit <- readRDS(r$fit_file); gk <- gwas_keys(fit)
  qp <- caqtl_prepare(copy(QC[[pk]]))
  o <- pair_test(gwas_row(fit, gk, r$idx1), gk, qp, fit$ld_rows, r$gwas_lead_snp)
  o <- o[!vapply(o, is.null, logical(1))]
  cbind(data.table(gwas = st, hit1 = h, peak = pk, refit_gene = r$refit_gene, idx1 = r$idx1, reproduced = r$reproduced,
                   gwas_lead_snp = r$gwas_lead_snp), as.data.table(o))
}
pc2 <- rbindlist(list(prim_ctrl("UKBB_GGT", "1:178514560", "peak21014"), prim_ctrl("UKBB_GGT", "5:166834153", "peak273749")), fill = TRUE)
cat("\n== PC2 (chr1, partly seen) and PC2b (chr5, blind) primary route ==\n"); print(pc2)
pass_prim <- function(x) isTRUE(x$testable) && isTRUE(x$PP.H4 > PP_ENTRY) && isTRUE(x$lead_ld_r2 >= LD_ENTRY)
S2a <- pc1[peak == "peak273749", PP.H4_p12_1e5] > 0.8
S2b <- pass_prim(pc2[1]); S2c <- pass_prim(pc2[2])
S2d <- abs(pc1[peak == "peak21014", PP.H4_p12_1e5] - 0.9466) <= 0.01 && abs(pc1[peak == "peak156441", PP.H4_p12_1e5] - 0.9990) <= 0.01
PC3 <- all(pc1$pc3_max_abs_diff < 1e-8)
S2_PASS <- isTRUE(S2a && S2b && S2c && S2d)
cat(sprintf("S2 (a) blind PC1 peak273749 PP.H4>0.8: %s | (b) PC2: %s | (c) PC2b: %s | (d) seen PC1 reproduce pilot: %s | PC3 (code check): %s\n",
            S2a, S2b, S2c, S2d, PC3))
cat("S2", if (S2_PASS) "PASS" else "FAIL (stop rule)", "\n")
fwrite(pc1, file.path(RES, "controls_pc1_pc3.tsv"), sep = "\t"); fwrite(pc2, file.path(RES, "controls_pc2_pc2b.tsv"), sep = "\t")
gate <- list(job = Sys.getenv("SLURM_JOB_ID"), S1_pass = S1_PASS, S1 = s1, S2_pass = S2_PASS,
             S2a = S2a, S2b = S2b, S2c = S2c, S2d = S2d, PC3 = PC3,
             selection_sha256 = readLines(file.path(DIAG, "selection_sha256.txt")))
jsonlite::write_json(gate, file.path(DIAG, "job1_gate.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
writeLines(capture.output(sessionInfo()), file.path(DIAG, "sessionInfo_job1.txt"))
stamp("H53 job 1 done. S1", S1_PASS, "S2", S2_PASS)
