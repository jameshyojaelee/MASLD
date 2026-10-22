# AN7 step 2: H59 primary (varIdent GLS, within-dataset stage permutation, BH over 117)
# rerun as two sensitivities. Uses the H59 library unchanged (sourced read-only from the
# H59 run dir): same fit_gls, perm_stage, key matrix (seed 20260923, 19,999 draws),
# donor universe (GSE202379 + GSE244832, 55 hepatocyte donors), G1 rule, perm p and
# selection rule (perm q < 0.05 over 117 and stability_median >= 0.5).
#   S1: donor scores = mean of stored per-cell scores over v2 QC-passing cells (job 1711 construction).
#   S2: published donor scores; T-cell donors with < 20 T cells removed; 11 T-cell programs refit;
#       BH over 117 with the other 106 programs' H59 perm p unchanged.
AN7 <- Sys.getenv("AN7_DIR")
SCORES <- Sys.getenv("AN7_SCORES", file.path(AN7, "out"))
H59 <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H59-20261015T112432Z"
Sys.setenv(H59_RUN_DIR = H59, H59_OUT = file.path(AN7, "out"))
source(file.path(H59, "scripts/h59_lib.R"))
stopifnot(OUT == file.path(AN7, "out"), dir.exists(OUT))
msg("AN7 step 2; B =", B, "cores =", NCORES, "OUT =", OUT)

# H59 check_hashes() also hashes the rendered Fig 4C PDF, which was redrawn after H59
# (display file, not an analysis input). Check every data/code input; record the PDF only.
check_hashes_an7 <- function() {
  ps <- sha256(file.path(RUN, "PRESPEC_H59.md"))
  if (ps != PRESPEC_SHA) stop("PRESPEC_H59.md sha256 mismatch: ", ps)
  obs <- vapply(INPUTS$role, function(r) sha256(inp(r)), "")
  res <- copy(INPUTS)[, observed := obs][, match := observed == sha256]
  if (!all(res[role != "fig4c_pdf", match]) || file.size(inp("scores_primary")) != 1869675) {
    print(res[match == FALSE]); stop("DATA INPUT HASH MISMATCH; refusing to run")
  }
  msg("data input hashes OK; fig4c_pdf match =", res[role == "fig4c_pdf", match])
  res
}
hashes <- check_hashes_an7()
wtsv(hashes, "inputs.sha256.tsv")
dat <- load_data()
programs <- dat$programs
pid <- paste(programs$cell_type, programs$module, sep = "_")
donors62 <- dat$donors62
K <- key_matrix(B)
msg("key matrix", nrow(K), "x", ncol(K), "first key", sprintf("%.10f", K[1, 1]))
prim <- fread(file.path(H59, "out/outcome/results_primary_family.tsv"))
setorder(prim, cell_type, module); stopifnot(identical(prim$program, pid))
Tn_h59 <- as.matrix(fread(file.path(H59, "out/null_gls_t.tsv.gz"))); stopifnot(identical(colnames(Tn_h59), pid))

perm_p <- function(obs, null) {   # verbatim from 02_h59_outcome.R
  vapply(seq_along(obs), function(j) {
    if (!is.finite(obs[j])) return(NA_real_)
    v <- abs(null[, j]); v <- v[is.finite(v)]
    (1 + sum(v >= abs(obs[j]) * (1 - 1e-12))) / (1 + length(v))
  }, 0)
}
g1_rule <- function(d) {          # G1 from 01_h59_gates_null.R, per program
  all(sapply(FIT_DATASETS, function(ds) { x <- d[dataset == ds, .N, by = stage]; sum(x$N) >= 8 && sum(x$N >= 2) >= 2 }))
}
attach_meta <- function(path) {
  x <- fread(path)
  x <- merge(x[, .(cell_type, module = as.integer(module), donor, score, n_cells)], dat$elig, by = "donor")
  x[dataset != DROP_DATASET & dataset %in% FIT_DATASETS]
}

# ------------------------------------------------ reproduction check (R0)
allc <- attach_meta(file.path(SCORES, "donor_scores_all_cells.tsv"))
pfit <- dat$p[dataset != DROP_DATASET & dataset %in% FIT_DATASETS]
chk <- merge(pfit[, .(cell_type, module, donor, score)], allc[, .(cell_type, module, donor, score_regen = score)],
             by = c("cell_type", "module", "donor"), all = TRUE)
r0_score <- max(abs(chk$score - chk$score_regen))
stopifnot(!anyNA(chk$score), !anyNA(chk$score_regen), r0_score < 1e-10)
r0 <- rbindlist(lapply(seq_len(nrow(programs)), function(i) {
  d <- prog_data(pfit, programs$cell_type[i], programs$module[i], FIT_DATASETS)
  f <- fit_gls(d$score, d$stage, d$dataset)
  data.table(program = pid[i], n_fit = nrow(d), gls_beta = f[["beta"]], gls_t = f[["t"]])
}))
r0[, perm_p := perm_p(gls_t, Tn_h59)][prim$g1_pass == FALSE, perm_p := 1][, perm_q := bh117(perm_p)]
r0chk <- data.table(check = c("max_abs_score_regen_vs_published", "max_abs_dbeta_vs_H59", "max_abs_dt_vs_H59",
                              "max_abs_dpermq_vs_H59", "max_abs_dn_vs_H59"),
                    value = c(r0_score, max(abs(r0$gls_beta - prim$gls_beta)), max(abs(r0$gls_t - prim$gls_t)),
                              max(abs(r0$perm_q - prim$perm_q)), max(abs(r0$n_fit - prim$n_fit))))
wtsv(r0chk, "R0_reproduce_H59.tsv")
print(r0chk)
stopifnot(r0chk$value[2] < 1e-8, r0chk$value[4] < 1e-12, r0chk$value[5] == 0)

# ------------------------------------------------ generic H59 family run
run_family <- function(scores, idx, label) {
  NL <- mclapply(idx, function(i) {
    d <- prog_data(scores, programs$cell_type[i], programs$module[i], FIT_DATASETS)
    obs <- fit_gls(d$score, d$stage, d$dataset)
    idp <- perm_stage(d, matrix(seq_len(62) / 63, nrow = 1), donors62)[1, ]
    P <- perm_stage(d, K, donors62)
    res <- t(vapply(seq_len(nrow(P)), function(k) fit_gls(d$score, P[k, ], d$dataset), numeric(6)))
    list(obs = obs, n = nrow(d), n202 = sum(d$dataset == "GSE202379"), n244 = sum(d$dataset == "GSE244832"),
         g1 = g1_rule(d), id_ok = identical(as.numeric(idp), as.numeric(d$stage)), null = res)
  }, mc.cores = NCORES, mc.preschedule = FALSE)
  bad <- which(vapply(NL, function(x) inherits(x, "try-error") || is.null(x), TRUE))
  if (length(bad)) stop(label, ": worker errors in ", paste(bad, collapse = ","))
  Tn <- sapply(NL, function(x) x$null[, "t"]); if (!is.matrix(Tn)) Tn <- matrix(Tn, nrow = 1)
  colnames(Tn) <- pid[idx]
  fwrite(as.data.table(Tn), file.path(OUT, paste0("null_gls_t_", label, ".tsv.gz")), sep = "\t")
  r <- data.table(program = pid[idx], cell_type = programs$cell_type[idx], module = programs$module[idx],
                  module_name = programs$module_name[idx], stability_median = programs$stability_median[idx],
                  n_fit = sapply(NL, `[[`, "n"), n_GSE202379 = sapply(NL, `[[`, "n202"), n_GSE244832 = sapply(NL, `[[`, "n244"),
                  g1_pass = sapply(NL, `[[`, "g1"), identity_stage_equal = sapply(NL, `[[`, "id_ok"),
                  gls_beta = sapply(NL, function(x) x$obs[["beta"]]), gls_se = sapply(NL, function(x) x$obs[["se"]]),
                  gls_t = sapply(NL, function(x) x$obs[["t"]]),
                  n_failed_null_draws = colSums(!is.finite(Tn)))
  r[, perm_p := perm_p(gls_t, Tn)]
  r[g1_pass == FALSE, perm_p := 1]
  r
}

# ------------------------------------------------ S1: QC-passing-cell donor scores, all 117
qcc <- attach_meta(file.path(SCORES, "donor_scores_v2qc_cells.tsv"))
msg("S1 start: QC rows", nrow(qcc))
t0 <- Sys.time()
s1 <- run_family(qcc, seq_len(nrow(programs)), "S1_v2qc")
msg("S1 null done in", format(Sys.time() - t0))
s1[, perm_q := bh117(perm_p)]
s1[, selected := is.finite(perm_q) & perm_q < 0.05 & is.finite(stability_median) & stability_median >= 0.5]
s1 <- merge(s1, prim[, .(program, h59_n_fit = n_fit, h59_gls_beta = gls_beta, h59_gls_t = gls_t,
                         h59_perm_p = perm_p, h59_perm_q = perm_q, h59_selected = selected_H59)], by = "program", sort = FALSE)
setorder(s1, cell_type, module)
wtsv(s1, "S1_v2qc_results_117.tsv")
ok <- is.finite(s1$gls_beta) & is.finite(s1$h59_gls_beta)
sp_all <- cor.test(s1$gls_beta[ok], s1$h59_gls_beta[ok], method = "spearman", exact = FALSE)
s1lin <- s1[ok, .(n_programs = .N, spearman_beta = cor(gls_beta, h59_gls_beta, method = "spearman"),
                  sign_flips = sum(sign(gls_beta) != sign(h59_gls_beta)),
                  median_abs_dbeta = median(abs(gls_beta - h59_gls_beta)), max_abs_dbeta = max(abs(gls_beta - h59_gls_beta)),
                  donors_lost_max = max(h59_n_fit - n_fit)), by = cell_type]
s1sum <- data.table(
  n_programs = nrow(s1), n_beta_finite = sum(ok), selected = sum(s1$selected), calls_q05 = sum(s1$perm_q < 0.05),
  smallest_perm_q = min(s1$perm_q), smallest_q_program = s1$program[which.min(s1$perm_q)],
  spearman_beta_vs_h59 = unname(sp_all$estimate), spearman_p = sp_all$p.value,
  spearman_t_vs_h59 = cor(s1$gls_t[ok], s1$h59_gls_t[ok], method = "spearman"),
  sign_flips = sum(sign(s1$gls_beta[ok]) != sign(s1$h59_gls_beta[ok])),
  median_abs_dbeta = median(abs(s1$gls_beta[ok] - s1$h59_gls_beta[ok])), max_abs_dbeta = max(abs(s1$gls_beta[ok] - s1$h59_gls_beta[ok])),
  g1_fail = sum(!s1$g1_pass), identity_fail = sum(!s1$identity_stage_equal),
  max_failed_draw_frac = max(s1$n_failed_null_draws) / nrow(K), programs_with_fewer_donors = sum(s1$n_fit < s1$h59_n_fit))
wtsv(s1sum, "S1_summary.tsv")
wtsv(s1lin, "S1_by_lineage.tsv")
print(s1sum); print(s1lin)
print(s1[program %in% c("hepatocytes_20", "hepatocytes_8")])

# ------------------------------------------------ S2: T-cell donors with < 20 T cells removed
TMIN <- 20L
tpub <- merge(pfit, allc[, .(cell_type, module, donor, n_cells)], by = c("cell_type", "module", "donor"))
stopifnot(nrow(tpub) == nrow(pfit))
gated <- tpub[!(cell_type == "tcells" & n_cells < TMIN)]
rm_t <- unique(tpub[cell_type == "tcells" & n_cells < TMIN, .(donor, dataset, stage, n_cells)])
wtsv(rm_t[order(dataset, n_cells)], "S2_tcell_donors_removed.tsv")
tidx <- which(programs$cell_type == "tcells")
stopifnot(length(tidx) == 11L)
msg("S2 start: T-cell donors removed", nrow(rm_t))
s2 <- run_family(gated, tidx, "S2_tcell_gate20")
fam <- copy(prim[, .(program, perm_p)])
fam[match(s2$program, program), perm_p := s2$perm_p]
fam[, perm_q := bh117(perm_p)]
s2[, perm_q_bh117 := fam$perm_q[match(program, fam$program)]]
s2[, perm_q_bh11 := p.adjust(perm_p, method = "BH")]
s2[, selected := is.finite(perm_q_bh117) & perm_q_bh117 < 0.05 & is.finite(stability_median) & stability_median >= 0.5]
s2 <- merge(s2, prim[, .(program, h59_n_fit = n_fit, h59_gls_beta = gls_beta, h59_perm_p = perm_p, h59_perm_q = perm_q)],
            by = "program", sort = FALSE)
setorder(s2, module)
wtsv(s2, "S2_tcell_gate20_results_11.tsv")
s2sum <- data.table(tcell_donors_removed = nrow(rm_t), removed_GSE202379 = sum(rm_t$dataset == "GSE202379"),
                    removed_GSE244832 = sum(rm_t$dataset == "GSE244832"), n_fit = paste(unique(s2$n_fit), collapse = ";"),
                    h59_n_fit = paste(unique(s2$h59_n_fit), collapse = ";"),
                    g1_fail = sum(!s2$g1_pass), selected = sum(s2$selected), calls_bh117_q05 = sum(s2$perm_q_bh117 < 0.05),
                    smallest_perm_p = min(s2$perm_p), smallest_perm_q_bh117 = min(s2$perm_q_bh117),
                    smallest_q_program = s2$program[which.min(s2$perm_q_bh117)], smallest_perm_q_bh11 = min(s2$perm_q_bh11),
                    family_calls_bh117_all = sum(fam$perm_q < 0.05),
                    h59_smallest_tcell_perm_p = min(prim[cell_type == "tcells", perm_p]),
                    h59_smallest_tcell_perm_q = min(prim[cell_type == "tcells", perm_q]))
wtsv(s2sum, "S2_summary.tsv")
print(s2sum); print(s2[, .(program, module_name, n_fit, g1_pass, gls_beta, h59_gls_beta, perm_p, perm_q_bh117, h59_perm_q)])

writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
file.create(file.path(OUT, "DONE"))
msg("AN7 step 2 complete")
