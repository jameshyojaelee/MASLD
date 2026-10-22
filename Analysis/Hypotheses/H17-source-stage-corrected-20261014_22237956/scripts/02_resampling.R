#!/usr/bin/env Rscript
# H17 step 2: label-permuted contrasts (negative control 2) and donor-level
# bootstrap of the Spearman correlation (S2).
# Usage: Rscript 02_resampling.R <HYP_DIR> <OUT_DIR>
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
hyp_dir <- args[1]; out_dir <- args[2]
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(here, "h17_lib.R"))
suppressPackageStartupMessages(library(parallel))
stopifnot(file.exists(file.path(out_dir, "STEP1_DONE")))
wt <- function(x, f) fwrite(x, file.path(out_dir, f), sep = "\t")
bind_ok <- function(l, what) {
  ok <- vapply(l, is.data.table, logical(1))
  if (!all(ok)) message(what, " failures: ", sum(!ok), "; first: ", as.character(l[!ok][[1]]))
  rbindlist(l[ok])
}
NC <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
N_PERM_B <- 1000L; N_SWAP_A <- 100L; N_BOOT <- 1000L

st <- readRDS(file.path(out_dir, "state_step1.rds"))
de_B <- fread(file.path(out_dir, "invivo_de.tsv"))
A <- load_ipsc(hyp_dir)
cA <- A$counts[st$genes_A, ]; shA <- A$sheet; rm(A); invisible(gc())
B <- load_invivo("steatotic")
cB <- B$counts[st$genes_B, ]; mB <- B$meta; rm(B); invisible(gc())
U <- st$U; idx_A <- st$idx_A; sA <- st$sA
qB_U <- de_B$q[match(U, de_B$gene)]; sB_U <- sign(de_B$logFC[match(U, de_B$gene)])

# ---------------------------------------------------------------- pre-generate indices
set.seed(H17_SEED)
perm_B <- lapply(seq_len(N_PERM_B), function(i) {
  g <- as.character(mB$group)
  for (ds in unique(mB$dataset)) { ix <- which(mB$dataset == ds); g[ix] <- g[ix][sample.int(length(ix))] }
  factor(g, levels = c("control", "steatotic"))
})
set.seed(H17_SEED)
parts <- levels(shA$participant)
swap_A <- lapply(seq_len(N_SWAP_A), function(i) parts[runif(length(parts)) < 0.5])
set.seed(H17_SEED)
boot_A <- lapply(seq_len(N_BOOT), function(i) sample(parts, length(parts), replace = TRUE))
strata <- split(seq_len(nrow(mB)), interaction(mB$dataset, mB$group, drop = TRUE))
boot_B <- lapply(seq_len(N_BOOT), function(i)
  unlist(lapply(strata, function(ix) ix[sample.int(length(ix), length(ix), replace = TRUE)]), use.names = FALSE))

# ---------------------------------------------------------------- in vivo label permutation
t0 <- Sys.time()
res_B <- mclapply(seq_len(N_PERM_B), function(i) {
  m <- copy(mB); m[, group := perm_B[[i]]]
  de <- fit_full(cB, st$genes_B, invivo_design(m), "groupsteatotic")
  lf <- de$logFC[match(U, de$gene)]; q <- de$q[match(U, de$gene)]
  sb <- sign(lf)
  data.table(perm = i, n_sig_invivo = sum(de$q < Q_THR),
             S1_agreement = mean(sA[idx_A] == sb[idx_A]),
             S2_spearman = cor(st$logFC_A, lf, method = "spearman"),
             n_sig_both = sum(st$q_A < Q_THR & q < Q_THR),
             n_core = sum(st$q_A < Q_THR & q < Q_THR & sA == sb))
}, mc.cores = NC, mc.preschedule = TRUE)
res_B <- bind_ok(res_B, "in vivo permutation")
wt(res_B, "negctrl_invivo_label_permutation.tsv")
message("in vivo label permutations: ", format(Sys.time() - t0))

# ---------------------------------------------------------------- iPSC within-pair label swaps
t0 <- Sys.time()
res_A <- mclapply(seq_len(N_SWAP_A), function(i) {
  tr <- as.character(shA$treatment)
  fl <- shA$participant %in% swap_A[[i]]
  tr[fl] <- ifelse(tr[fl] == "HepC", "HepFA", "HepC")
  tr <- factor(tr, levels = c("HepC", "HepFA"))
  de <- fit_full(cA, st$genes_A, ipsc_design(shA$participant, tr), "treatmentHepFA")
  lf <- de$logFC[match(U, de$gene)]; q <- de$q[match(U, de$gene)]
  data.table(swap = i, n_swapped = length(swap_A[[i]]), n_sig_ipsc = sum(de$q < Q_THR),
             n_core = sum(q < Q_THR & qB_U < Q_THR & sign(lf) == sB_U))
}, mc.cores = NC, mc.preschedule = TRUE)
res_A <- bind_ok(res_A, "iPSC swap")
wt(res_A, "negctrl_ipsc_label_swap.tsv")
message("iPSC label swaps: ", format(Sys.time() - t0))

# ---------------------------------------------------------------- bootstrap S2
t0 <- Sys.time()
lib_of <- split(seq_len(nrow(shA)), shA$participant)
res_boot <- mclapply(seq_len(N_BOOT), function(b) {
  ps <- boot_A[[b]]
  cols <- unlist(lapply(ps, function(p) lib_of[[p]]), use.names = FALSE)
  pp <- factor(rep(seq_along(ps), times = vapply(ps, function(p) length(lib_of[[p]]), integer(1))))
  ca <- cA[, cols]; colnames(ca) <- make.unique(colnames(ca))
  lfA <- fit_logfc(ca, st$genes_A, ipsc_design(pp, shA$treatment[cols]), "treatmentHepFA")
  ib <- boot_B[[b]]
  cb <- cB[, ib]; colnames(cb) <- make.unique(colnames(cb))
  lfB <- fit_logfc(cb, st$genes_B, invivo_design(mB[ib]), "groupsteatotic")
  data.table(boot = b, S2_spearman = cor(lfA[U], lfB[U], method = "spearman"))
}, mc.cores = NC, mc.preschedule = TRUE)
res_boot <- bind_ok(res_boot, "bootstrap")
wt(res_boot, "bootstrap_S2.tsv")
message("bootstrap: ", format(Sys.time() - t0))
writeLines("ok", file.path(out_dir, "STEP2_DONE"))
