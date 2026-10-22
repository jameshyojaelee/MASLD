#!/usr/bin/env Rscript
# H36 step 5: calls, permutation calibration, direction consistency, composition
# sensitivity (S1), leave-one-cohort-out (descriptive), and concordance with sQTL
# colocalization. PRESPEC_H36.md sections 4, 5, 7 and Addenda 1 to 3.
# Usage: Rscript h36_05_integrate.R <repo> <h36_dir> <out_dir> [stage_dir_name=stage_test] [force_calls=0]
# force_calls > 0 is a code-path test switch: the n clusters with the smallest p are treated as calls (never used for results).
suppressPackageStartupMessages({ library(data.table); library(parallel); library(GenomicRanges); library(rtracklayer) })
args <- commandArgs(trailingOnly = TRUE)
REPO <- args[1]; H36 <- args[2]; OUT <- args[3]
STAGE_DIR <- if (length(args) >= 4) args[4] else "stage_test"
FORCE_CALLS <- if (length(args) >= 5) as.integer(args[5]) else 0L
NC <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))
HERE <- file.path(REPO, "scripts/analysis/hypotheses_2026q4/H36_splicing")
source(file.path(HERE, "h36_dm.R")); source(file.path(HERE, "h36_prep.R"))
set.seed(20260923)
BUILD <- file.path(H36, "build"); ST <- file.path(H36, STAGE_DIR); SQ <- file.path(H36, "sqtl_coloc")
GTF <- "/gpfs/commons/home/jameslee/reference_genome/cellranger/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz"
COMP <- file.path(REPO, "Analysis/Deconvolution/corrected_774_20261009T141331Z/summary/composition_new_774.tsv.gz")
load_build(BUILD)
w <- function(x, f) fwrite(x, file.path(OUT, f), sep = "\t")

# ---- 1. observed calls ----------------------------------------------------------
obs <- fread(file.path(ST, "stage_test_observed.tsv.gz"))
perm <- fread(file.path(ST, "stage_test_perm.tsv.gz"))
status_tab <- obs[, .N, by = .(status = sub(":.*", "", status))]
w(status_tab, "cluster_status_counts.tsv"); print(status_tab)
ok <- obs[status == "ok"]
perm <- perm[match(ok$cl, cl)]
stopifnot(identical(perm$cl, ok$cl))
K <- sum(grepl("^lr_perm", names(perm)))
ok[, q := p.adjust(p, "BH")]
ok[, conv_ok := conv0 == 0 & conv1 == 0]
ok[, call_noconv := q < 0.05 & max_abs_dpsi >= 0.10]
ok[, call_bh := call_noconv & conv_ok]
cat(sprintf("BH family: %d clusters; q<0.05: %d; q<0.05 and max|dPSI|>=0.10: %d (%d with both fits converged); non-converged fits: %d\n",
            nrow(ok), sum(ok$q < 0.05), sum(ok$call_noconv), sum(ok$call_bh), sum(!ok$conv_ok)))

# ---- 2. permutation calibration and the call rule fixed in PRESPEC section 5 ------------
P <- sapply(seq_len(K), function(k) pchisq(perm[[paste0("lr_perm", k)]], ok$df, lower.tail = FALSE))
D <- sapply(seq_len(K), function(k) perm[[paste0("dpsi_perm", k)]])
CV <- sapply(seq_len(K), function(k) perm[[paste0("conv_perm", k)]] == 0 & ok$conv0 == 0)
calls_perm_bh <- lapply(seq_len(K), function(k) which(p.adjust(P[, k], "BH") < 0.05 & D[, k] >= 0.10 & CV[, k]))
cal <- data.table(set = c("observed", paste0("perm", seq_len(K))),
                  n_clusters = nrow(ok),
                  frac_p_lt_0.05 = c(mean(ok$p < 0.05), colMeans(P < 0.05)),
                  n_p_lt_1e4 = c(sum(ok$p < 1e-4), colSums(P < 1e-4)),
                  n_q_lt_0.05 = c(sum(ok$q < 0.05), apply(P, 2, function(x) sum(p.adjust(x, "BH") < 0.05))),
                  n_calls_bh = c(sum(ok$call_bh), lengths(calls_perm_bh)))
cal[, `:=`(n_converged = c(sum(ok$conv_ok), colSums(CV)),
            frac_p_lt_0.05_converged = c(mean(ok$p[ok$conv_ok] < 0.05), vapply(seq_len(K), function(k) mean(P[CV[, k], k] < 0.05), numeric(1))),
            frac_p_lt_0.05_not_converged = c(mean(ok$p[!ok$conv_ok] < 0.05), vapply(seq_len(K), function(k) mean(P[!CV[, k], k] < 0.05), numeric(1))))]
w(cal, "calibration_by_set.tsv"); print(cal)
mean_frac <- mean(cal[set != "observed"]$frac_p_lt_0.05)
bh_stands <- mean_frac >= 0.03 && mean_frac <= 0.07
cat(sprintf("mean permuted fraction p<0.05 = %.4f; accepted range 0.03-0.07 -> %s\n", mean_frac,
            if (bh_stands) "BH rule stands" else "empirical-FDR rule replaces it"))
# empirical FDR rule: FDR(t) = mean_k #{p_perm_k <= t} / #{p_obs <= t}; calls at the largest t with FDR < 0.05
ps <- sort(ok$p); allperm <- sort(as.vector(P))
fdr_t <- (findInterval(ps, allperm) / K) / seq_along(ps)
tstar <- if (any(fdr_t < 0.05)) max(ps[fdr_t < 0.05]) else NA_real_
ok[, call_emp := !is.na(tstar) & p <= tstar & max_abs_dpsi >= 0.10]
cat(sprintf("empirical-FDR rule: t* = %s; clusters p<=t* with |dPSI|>=0.10: %d\n", format(tstar, digits = 3), sum(ok$call_emp)))
ok[, call := if (bh_stands) call_bh else call_emp]
if (FORCE_CALLS > 0L) { cat("CODE-PATH TEST MODE: forcing", FORCE_CALLS, "calls\n"); ok[order(p)[seq_len(FORCE_CALLS)], call := TRUE] }
calls_perm <- if (bh_stands) calls_perm_bh else lapply(seq_len(K), function(k) which(P[, k] <= tstar & D[, k] >= 0.10 & CV[, k]))
cat(sprintf("PRIMARY CALLS (%s): %d of %d clusters\n", if (bh_stands) "BH" else "empirical FDR", sum(ok$call), nrow(ok)))

# ---- 3. annotate clusters with gene symbols (cluster span, same strand) ----------------
gtf <- import(GTF, feature.type = "gene")
ann_genes <- function(chr, s, e, str) {
  g <- GRanges(chr, IRanges(s, e), strand = str)
  ov <- findOverlaps(g, gtf, ignore.strand = FALSE)
  tapply(gtf$gene_name[subjectHits(ov)], queryHits(ov), function(x) paste(head(unique(x), 3), collapse = ","))
}
cn <- ok[, .(cl, chr, str, start_min, end_max)]
cn[, strand := ifelse(str %in% c("+", "-"), str, "*")]
gn <- ann_genes(cn$chr, cn$start_min, cn$end_max, cn$strand)
ok[, genes := NA_character_]; ok[as.integer(names(gn)), genes := unname(gn)]

# ---- 4. direction consistency, S1 composition, leave-one-cohort-out (call clusters) --------
comp <- fread(COMP)[, .(sample_id, hep = Hepatocytes, mac = Macrophages, fib = Fibroblasts)]
pc <- merge(pers[, .(pi, sample_id)], comp, by = "sample_id", all.x = TRUE)[order(pi)]
stopifnot(identical(pc$pi, pers$pi), !anyNA(pc$hep))
Z <- scale(as.matrix(pc[, .(hep, mac, fib)]))     # z-scores over all 774 persons
sgn <- function(x) sign(x)
extra_checks <- function(k) {
  r <- ok[cl == k]
  p <- prep(k, TRUE); Y <- p$Y; pi <- p$pi
  st <- stage0[pi]; d <- ds0[pi]
  col <- which(p$jr == r$top_jr); y <- Y[, col] / rowSums(Y)
  rho <- vapply(sort(unique(d)), function(dd) {
    i <- d == dd
    if (sum(i) < 30L || sd(st[i]) == 0 || sd(y[i]) == 0) NA_real_ else cor(st[i], y[i], method = "spearman")
  }, numeric(1))
  ev <- !is.na(rho); same <- sum(sgn(rho[ev]) == sgn(r$top_dpsi))
  # S1: add composition z-scores to the null and stage models
  X0 <- design0(pi); X0c <- independent_cols(cbind(X0, Z[pi, , drop = FALSE]))
  xs <- st - mean(st)
  s1 <- dm_lrt(Y, X0c, xs); xb <- colMeans(X0c); m <- mean(st)
  dp1 <- as.numeric(dm_pred(s1$f1, c(xb, 4 - m)) - dm_pred(s1$f1, c(xb, 0 - m)))[col]
  # leave one cohort out: nominal p < 0.05 and same sign of the top-junction dPSI
  loco <- vapply(sort(unique(d)), function(dd) {
    i <- d != dd
    if (length(unique(d[i])) < 2L || sum(i) < 100L) return(NA_character_)
    Xl <- independent_cols(design0(pi[i])); xl <- st[i] - mean(st[i])
    rl <- dm_lrt(Y[i, , drop = FALSE], Xl, xl); xbl <- colMeans(Xl); ml <- mean(st[i])
    dpl <- as.numeric(dm_pred(rl$f1, c(xbl, 4 - ml)) - dm_pred(rl$f1, c(xbl, 0 - ml)))[col]
    paste0(signif(rl$p, 3), "|", signif(dpl, 3))
  }, character(1))
  loco_p <- as.numeric(sub("\\|.*", "", loco)); loco_d <- as.numeric(sub(".*\\|", "", loco))
  data.table(cl = k, n_evaluable_cohorts = sum(ev), n_same_sign = same,
             consistent = sum(ev) >= 3L && same == sum(ev),
             cohort_rho = paste(names(rho), signif(rho, 2), sep = ":", collapse = ";"),
             s1_p = s1$p, s1_top_dpsi = dp1, s1_keep = s1$p < 0.05 && sgn(dp1) == sgn(r$top_dpsi),
             loco_each = paste(names(loco), loco, sep = "=", collapse = ";"),
             loco_n_eval = sum(!is.na(loco_p)),
             loco_n_keep = sum(loco_p < 0.05 & sgn(loco_d) == sgn(r$top_dpsi), na.rm = TRUE))
}
cl_calls <- ok[call == TRUE, cl]
chk <- if (length(cl_calls)) rbindlist(mclapply(cl_calls, function(k) tryCatch(extra_checks(k), error = function(e) data.table(cl = k, consistent = NA)),
                          mc.cores = NC, mc.preschedule = TRUE), fill = TRUE) else
  data.table(cl = integer(), n_evaluable_cohorts = integer(), n_same_sign = integer(), consistent = logical(), s1_p = numeric(),
             s1_top_dpsi = numeric(), s1_keep = logical(), loco_n_eval = integer(), loco_n_keep = integer())
ok <- merge(ok, chk, by = "cl", all.x = TRUE)
setorder(ok, p)
w(ok[, !c("dpsi_jr", "cohorts")][, .(cl, chr, str, start_min, end_max, genes, J_used, n_persons, lr, df, p, q, max_abs_dpsi, top_jr, top_dpsi,
                                    conv0, conv1, call_bh, call_emp, call, n_evaluable_cohorts, n_same_sign, consistent, s1_p, s1_top_dpsi, s1_keep,
                                    loco_n_eval, loco_n_keep)], "clusters_tested.tsv.gz")
w(ok[call == TRUE][order(q, -max_abs_dpsi)], "calls_primary_full.tsv")
cat(sprintf("calls: %d | direction-consistent in >=3 cohorts: %d of %d evaluable | S1 keep: %d of %d | LOCO keep in all evaluable leave-outs: %d\n",
            sum(ok$call), sum(ok$consistent, na.rm = TRUE), sum(!is.na(ok$consistent) & ok$call), sum(ok$s1_keep, na.rm = TRUE),
            sum(!is.na(ok$s1_keep)), sum(ok$call & ok$loco_n_eval > 0 & ok$loco_n_keep == ok$loco_n_eval, na.rm = TRUE)))
print(head(ok[call == TRUE, .(cl, chr, genes, J_used, n_persons, p = signif(p, 3), q = signif(q, 3), top_dpsi = signif(top_dpsi, 3), consistent, s1_keep)], 25))

# ---- 5. pipeline check: cohort effect --------------------------------------------------
cc <- fread(file.path(ST, "cohort_effect_pipeline_check.tsv"))
cat(sprintf("cohort-effect pipeline check: %d clusters, median p %.3g, fraction p<1e-4 %.3f\n", nrow(cc), median(cc$p), mean(cc$p < 1e-4)))

# ---- 6. concordance with colocalized sQTL introns ----------------------------------------
gate_map <- readLines(file.path(SQ, "GATE_mapping.txt"))[1]
cq <- fread(file.path(SQ, "sqtl_coloc_locus_intron.tsv.gz"))
cq[, tested_cl := !is.na(cl) & cl %in% ok$cl]
cq[, coloc := PP.H4 > 0.5]
n_coloc_introns_mapped <- uniqueN(cq[coloc == TRUE & tested_cl == TRUE]$molecular_trait_id)
stop_rule <- gate_map != "PASS" || n_coloc_introns_mapped < 20L
cat(sprintf("mapping gate: %s | colocalized introns mapping to a tested cluster: %d | stop rule (mapping fail or <20): %s\n",
            gate_map, n_coloc_introns_mapped, stop_rule))
pairs <- cq[tested_cl == TRUE]
conc_loci <- function(call_cl) unique(pairs[coloc == TRUE & cl %in% call_cl]$locus_id)
obs_conc <- conc_loci(ok[call == TRUE]$cl)
perm_conc <- vapply(calls_perm, function(ix) length(conc_loci(ok$cl[ix])), integer(1))
conc <- data.table(quantity = c("loci with a tested, colocalized sQTL intron mapped to a tested cluster",
                                "concordant loci (observed)", "concordant loci (permutation mean)",
                                "concordant loci (permutation max)", "permutation sets >= observed", "empirical p = (1 + n_ge) / (K + 1)"),
                   value = c(uniqueN(pairs[coloc == TRUE]$locus_id), length(obs_conc), mean(perm_conc), max(perm_conc),
                             sum(perm_conc >= length(obs_conc)), (1 + sum(perm_conc >= length(obs_conc))) / (K + 1)))
w(conc, "concordance_summary.tsv"); print(conc)
w(pairs[coloc == TRUE & cl %in% ok[call == TRUE]$cl][order(-PP.H4)][, .(locus_id, liver_study, molecular_trait_id, gene_id, cl, PP.H4, min_qtl_p, nsnps)], "concordant_locus_intron_pairs.tsv")
sym <- setNames(gtf$gene_name, sub("\\..*$", "", gtf$gene_id))
col_tab <- cq[coloc == TRUE][order(-PP.H4)]
col_tab <- merge(col_tab, ok[, .(cl, cl_genes = genes, stage_p = p, stage_q = q, max_abs_dpsi, top_dpsi, stage_call = call)], by = "cl", all.x = TRUE, sort = FALSE)
col_tab[, sqtl_gene := unname(sym[sub("\\..*$", "", gene_id)])]
w(col_tab[order(-PP.H4), .(locus_id, liver_study, molecular_trait_id, sqtl_gene, PP.H4, PP.H3, min_qtl_p, nsnps, map_class, cl, cl_genes,
                           stage_p, stage_q, max_abs_dpsi, top_dpsi, stage_call)], "colocalized_introns_stage_status.tsv")
# odds ratio: colocalized versus not, in a stage-associated cluster versus not; unit = (locus, intron) pair, bootstrap over loci
or_counts <- function(call_cl) {
  d <- pairs[, .(locus_id, coloc, call = cl %in% call_cl)]
  d[, .(a = sum(coloc & call), b = sum(coloc & !call), c = sum(!coloc & call), d = sum(!coloc & !call)), by = locus_id]
}
or_est <- function(m) { a <- m[1] + 0.5; b <- m[2] + 0.5; c <- m[3] + 0.5; d <- m[4] + 0.5; (a * d) / (b * c) }
M <- or_counts(ok[call == TRUE]$cl); Mm <- as.matrix(M[, .(a, b, c, d)])
boot <- replicate(1000L, or_est(colSums(Mm[sample.int(nrow(Mm), replace = TRUE), , drop = FALSE])))
perm_or <- vapply(calls_perm, function(ix) if (!length(ix)) NA_real_ else or_est(colSums(as.matrix(or_counts(ok$cl[ix])[, .(a, b, c, d)]))), numeric(1))
orr <- data.table(n_pairs = nrow(pairs), n_loci = uniqueN(pairs$locus_id), a_coloc_call = sum(Mm[, 1]), b_coloc_nocall = sum(Mm[, 2]),
                  c_nocoloc_call = sum(Mm[, 3]), d_nocoloc_nocall = sum(Mm[, 4]), OR_haldane = or_est(colSums(Mm)),
                  ci_lo = quantile(boot, 0.025), ci_hi = quantile(boot, 0.975), perm_OR_mean = if (all(is.na(perm_or))) NA_real_ else mean(perm_or, na.rm = TRUE),
                  perm_OR_max = if (all(is.na(perm_or))) NA_real_ else max(perm_or, na.rm = TRUE))
# the odds ratio is undefined without pairs in a stage-associated cluster or without colocalized pairs
if (orr$a_coloc_call + orr$c_nocoloc_call == 0 || orr$a_coloc_call + orr$b_coloc_nocall == 0) orr[, `:=`(OR_haldane = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_)]
w(orr, "odds_ratio_summary.tsv"); print(t(orr))
writeLines(c(sprintf("stop_rule_triggered=%s", stop_rule), sprintf("mapping_gate=%s", gate_map),
             sprintf("n_coloc_introns_mapped_to_tested_cluster=%d", n_coloc_introns_mapped)), file.path(OUT, "stop_rule.txt"))
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo_05.txt"))
cat("done\n")
