# Debug only: fabricate SYNTHETIC job 1-3 outputs with the exact schemas of 01/02/03 (states drawn from a
# made-up logistic curve, not from any fit), so 04_assemble.R can be exercised end to end before the real
# outcome exists. Nothing here is a result. Usage: Rscript debug_fake_outputs.R <debug_setup_dir> <out_root>
H80 <- "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H80-20261015T130802Z"
Sys.setenv(H80_SCRIPTS = file.path(H80, "scripts"))
source(file.path(H80, "scripts/h80_lib.R"))
a <- commandArgs(trailingOnly = TRUE); SET <- a[1]; ROOTD <- a[2]
fs <- file.path(ROOTD, "fake_setup"); fa <- file.path(ROOTD, "fake_armA"); fb <- file.path(ROOTD, "fake_armB_pc1"); fr <- file.path(ROOTD, "fake_readout")
for (d in c(fs, fa, fb, fr)) dir.create(d)
for (f in c("windows.tsv", "pc1_windows.tsv", "streams.rds", "setup_status.tsv", "pc2_results.tsv", "pc3_results.tsv",
            "D4_notconverged_summary_genes.tsv", "gene_state_by_fit.tsv", "input_hash_check.tsv")) file.copy(file.path(SET, f), file.path(fs, f))
writeLines("PASS", file.path(fs, "STATUS"))
W <- fread(file.path(fs, "windows.tsv")); P1 <- fread(file.path(fs, "pc1_windows.tsv"))
set.seed(1)
ppos <- function(le, lg) plogis(1.2 * (pmin(le, lg) - 5))
states_for <- function(n, p) {
  u <- runif(n); s <- ifelse(u < p, "positive", sample(c("no_eqtl_cs", "tested_h3", "tested_other", "untestable"), n, TRUE, c(.5, .3, .15, .05)))
  pp4 <- ifelse(s == "positive", runif(n, .5, 1), ifelse(s %in% c("tested_h3", "tested_other"), runif(n, 0, .5), NA))
  pp3 <- ifelse(s == "tested_h3", runif(n, .8, 1), ifelse(s %in% c("positive", "tested_other"), runif(n, 0, .3), NA))
  data.table(state = s, pp4 = pp4, pp3 = pp3, n_pairs = 1L, n_eligible = 1L, elig_match = TRUE, coloc_drew = FALSE, coloc_restored = TRUE)
}
cells <- list(); fits <- list(); reps <- list(); dg <- list()
cl <- rbind(data.table(cell_set = "E", lam_e = A_ESWEEP, lam_g = 10, causal = "shared"),
            data.table(cell_set = "G", lam_e = 10, lam_g = A_GSWEEP, causal = "shared"),
            data.table(cell_set = "D", lam_e = c(6, 10, 20), lam_g = c(10, 10, 20), causal = "distinct"))
for (r in W$row) for (rep in 1:10) {
  na_k <- r %% 37 == 0
  x <- cl[, .(row = r, rep = rep, cell_set, lam_e, lam_g, causal)]
  p <- ifelse(x$causal == "shared", ppos(x$lam_e, x$lam_g), 0.02)
  x <- cbind(x, states_for(nrow(x), p))
  if (na_k) x[causal == "distinct", state := "distinct_not_applicable"]
  cells[[length(cells) + 1L]] <- x
  f <- rbind(data.table(side = "eqtl", lambda = A_EQTL_LAMBDAS, causal_type = "shared"),
             data.table(side = "gwas", lambda = A_GWAS_SHARED, causal_type = "shared"),
             if (!na_k) data.table(side = "gwas", lambda = A_GWAS_DISTINCT, causal_type = "distinct"))
  f[, `:=`(row = r, rep = rep, error = NA_character_, converged = TRUE, n_cs = rbinom(.N, 1, plogis(1.2 * (lambda - 5))), cs_sizes = "",
           niter = 10L, lead_absz = 3.5 + 0.9 * lambda + rnorm(.N, 0, .5), z_causal = lambda + rnorm(.N), drew_rng = runif(.N) < .3,
           restored = TRUE, time_s = runif(.N, 5, 15))]
  fits[[length(fits) + 1L]] <- f
  reps[[length(reps) + 1L]] <- data.table(row = r, rep = rep, c_id = paste0("1:", sample(1e6, 1)), k_id = if (na_k) NA_character_ else "1:5",
                                          pool_size = if (na_k) 5L else 900L, epsE_hash = digest(runif(1)), epsG_hash = digest(runif(1)),
                                          noise_at_c = rnorm(1), stream_unchanged_by_fits = TRUE, kind_after = "L'Ecuyer-CMRG")
}
for (r in W$row) dg[[length(dg) + 1L]] <- data.table(row = r, ensembl = W[row == r, ensembl], stratum = W[row == r, stratum], ancestry = W[row == r, ancestry],
                                                     pE = W[row == r, nE], pG = 100L, pEG = 100L, ids_E_identical_rds = TRUE, n_dup_ids_E = 0L, na_E = 0L, na_G = 0L,
                                                     min_eig_E = 1e-3, n_clip_E = 0L, n_blocks_E = 1L, min_eig_G = 1e-3, n_clip_G = 0L, n_blocks_G = 1L,
                                                     trS_E = 1, trS2_E = 1, eff_df_E = 500, T_w = rnorm(1), mean_noise_c = rnorm(1, 0, .3),
                                                     nc3_T_pass = TRUE, nc3_causal_pass = TRUE, t_load_s = 1, t_eig_s = 1, t_total_s = 100, peak_rss_gb = 2, files_used = "x")
fwrite(rbindlist(cells), file.path(fa, "cells.tsv.gz"), sep = "\t"); fwrite(rbindlist(fits), file.path(fa, "fits.tsv.gz"), sep = "\t")
fwrite(rbindlist(reps), file.path(fa, "reps.tsv.gz"), sep = "\t"); fwrite(rbindlist(dg), file.path(fa, "window_diagnostics.tsv"), sep = "\t")
fwrite(data.table(row = W$row, ok = TRUE, error = NA_character_, t_total_s = 100), file.path(fa, "window_status.tsv"), sep = "\t")
writeLines("DONE", file.path(fa, "DONE"))
# Arm B
cb <- list(); fb_ <- list(); rb <- list(); db <- list()
for (r in W[armB == TRUE, row]) {
  for (rep in 1:10) for (gen in c("matched", "mismatched")) {
    x <- rbind(data.table(cell_set = "E", lam_e = B_ESWEEP, lam_g = 10, causal = "shared"), data.table(cell_set = "D", lam_e = 10, lam_g = 10, causal = "distinct"))
    p <- ifelse(x$causal == "shared", ppos(x$lam_e, 10) - (gen == "mismatched") * 0.05, 0.03)
    cb[[length(cb) + 1L]] <- cbind(data.table(row = r, rep = rep, generator = gen), x, states_for(nrow(x), pmax(p, 0)))
    fb_[[length(fb_) + 1L]] <- data.table(row = r, rep = rep, generator = gen, side = c(rep("eqtl", 6), "gwas", "gwas"), lambda = c(B_EQTL_LAMBDAS, 10, 10),
                                          causal_type = c(rep("shared", 7), "distinct"), error = NA_character_, converged = TRUE, n_cs = 1L, cs_sizes = "",
                                          niter = 5L, lead_absz = 5, z_causal = 1, drew_rng = FALSE, restored = TRUE, time_s = 5)
  }
  for (rep in 1:10) rb[[length(rb) + 1L]] <- data.table(row = r, rep = rep, c_id = paste0("1:", rep), k_id = "1:9", pool_size = 500L,
                                                        epsE_hash = digest(runif(1)), epsG_hash = "x", stream_unchanged_by_fits = TRUE, kind_after = "L'Ecuyer-CMRG")
  db[[length(db) + 1L]] <- data.table(arm = "B", row = r, T_w_matched = rnorm(1), T_w_mismatched = rnorm(1), files_used = "y")
}
fwrite(rbindlist(cb), file.path(fb, "armB_cells.tsv.gz"), sep = "\t"); fwrite(rbindlist(fb_), file.path(fb, "armB_fits.tsv.gz"), sep = "\t")
fwrite(rbindlist(rb), file.path(fb, "armB_reps.tsv.gz"), sep = "\t"); fwrite(rbindlist(db), file.path(fb, "armB_window_diagnostics.tsv"), sep = "\t")
cp <- list(); fp <- list(); rp <- list()
for (r in P1$row) for (rep in 1:10) {
  cp[[length(cp) + 1L]] <- cbind(data.table(row = r, rep = rep, lam_e = 12, lam_g = 20), states_for(1, 0.95))
  fp[[length(fp) + 1L]] <- data.table(row = r, rep = rep, side = c("eqtl", "gwas"), lambda = c(12, 20), causal_type = "shared", error = NA_character_,
                                      converged = TRUE, n_cs = 1L, cs_sizes = "", niter = 5L, lead_absz = 12, z_causal = 12, drew_rng = FALSE, restored = TRUE, time_s = 5)
  rp[[length(rp) + 1L]] <- data.table(row = r, rep = rep, epsE_hash = digest(runif(1)), stream_unchanged_by_fits = TRUE, kind_after = "L'Ecuyer-CMRG")
}
fwrite(rbindlist(cp), file.path(fb, "pc1_cells.tsv.gz"), sep = "\t"); fwrite(rbindlist(fp), file.path(fb, "pc1_fits.tsv.gz"), sep = "\t")
fwrite(rbindlist(rp), file.path(fb, "pc1_reps.tsv.gz"), sep = "\t")
fwrite(data.table(task = c(paste0("B", W[armB == TRUE, row]), paste0("P", P1$row)), ok = TRUE, error = NA_character_, t_total_s = 1), file.path(fb, "window_status.tsv"), sep = "\t")
writeLines("DONE", file.path(fb, "DONE"))
# readout
gx <- fread(file.path(fs, "gene_state_by_fit.tsv"))
L <- gx[fit %in% c("cs0", "cs_ge1"), .(ensembl, chr = chr_fit, rds = "x", nE = pmax(n_snps, 50L), rds_converged = TRUE, rds_n_cs = as.integer(fit == "cs_ge1"),
                                      lead_absz = ifelse(fit == "cs_ge1", rnorm(.N, as.numeric(Sys.getenv("FAKE_CSG_MEAN", "9")), 2), rnorm(.N, 4, 1.5)), n_pos_matched = n_snps, n_dup_pos = rbinom(.N, 1, .05),
                                      susie_state_t12, fit)]
L[, lead_absz_nodup := lead_absz]
fwrite(L, file.path(fr, "lead_z.tsv"), sep = "\t"); writeLines("DONE", file.path(fr, "DONE"))
cat("fake outputs written under", ROOTD, "\n")
